library(dplyr)
library(ggplot2)
library(cowplot)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
theme_update(panel.grid.major=element_line(size=0.5, color="grey", linetype="dashed"))

#####################################################
## Simulator of time-to-event data for any empirical hazard function h(t)
## --  clean version  --
## Slow, but does not require to know the inverse cumulative hazard fn.
## Currently simulates pairs of full siblings, to investigate changes in covariance
## after introducing mutations with different effects, timing and so on.
##                                              --- Julius, 2016-04-18

snpTableToMatrix = function(SnpTable){
        ## converts a table with numbers of SNPs of different classes
        ## and output a matrix with the effects of individual SNPs
        SnpMatrix = data.frame(name = rep(SnpTable$name, SnpTable$number),
                               effect = rep(SnpTable$effect, SnpTable$number),
                               spread = rep(SnpTable$spread, SnpTable$number),
                               maf = rep(SnpTable$maf, SnpTable$number),
                               time = rep(SnpTable$time, SnpTable$number),
                               type = rep(SnpTable$type, SnpTable$number))
        
        return(SnpMatrix)
}

snpMatrixToBetas = function(SnpMatrix, maxTime){
        ## generates beta coefficients (effect sizes) for each genomic position and time point

        BetaMatrix = matrix(0, nrow=nrow(SnpMatrix), ncol=maxTime-150)
        for(s in 1:nrow(SnpMatrix)){
                type = SnpMatrix$type[s]
                b = SnpMatrix$effect[s]
                
                if(type=="normal"){
                        m = SnpMatrix$time[s]
                        sd = SnpMatrix$spread[s]
                        BetaMatrix[s,] = dnorm(1:(maxTime-150), mean=m, sd=sd)*b
                } else if(type=="uniform"){
                        BetaMatrix[s,] = rep(b, maxTime-150)
                } else if(type=="stepon"){
                        t = SnpMatrix$time[s]
                        BetaMatrix[s,] = c(rep(0, t-1), rep(b, maxTime-150-t+1))
                } else if(type=="stepoff"){
                        t = SnpMatrix$time[s]
                        BetaMatrix[s,] = c(rep(b, t-1), rep(0, maxTime-150-t+1))
                }
        }
        ## remove negligible betas to speed up things and save memory
        BetaMatrix[abs(BetaMatrix)<1e-3] = 0
        BetaMatrix = Matrix(BetaMatrix)
        return(BetaMatrix)        
}

cppFunction("arma::mat generateGenome(int inds, arma::vec maf){
        arma::mat out = arma::randu(inds, maf.size());
            
        for(int g=0; g<maf.size(); g++){
            for(int i=0; i<inds; i++){
                   if(out(i, g) < maf(g)*maf(g)){
                        out(i, g) = 2;
                   } else if(out(i, g) < maf(g)*(2-maf(g))){
                        out(i, g) = 1;
                   } else {
                        out(i, g) = 0;
                   }
                }
            }
       return(out);
}", depends="RcppArmadillo")

cppFunction("arma::mat generateSibGenome(arma::mat genome, arma::vec maf){
        arma::mat out = arma::randu(genome.n_rows, maf.size());
            
        for(int g=0; g<maf.size(); g++){
            for(int i=0; i<genome.n_rows; i++){
                if(genome(i, g) == 0){
                        if( out(i, g) < maf(g)*maf(g)/4 ){
                                out(i, g) = 2;
                        } else if( out(i, g) < maf(g)*maf(g)*3/4 + maf(g)*(1-maf(g)) ){
                                out(i, g) = 1;
                        } else {
                                out(i, g) = 0;
                        }
                } else if(genome(i, g) == 1){
                        if( out(i, g) < maf(g)*(1-maf(g))/4 + maf(g)*maf(g)/2 ){
                                out(i, g) = 2;
                        } else if( out(i, g) < (1 - maf(g)*(1-maf(g)))/2 ){
                                out(i, g) = 0;
                        } else {
                                out(i, g) = 1;
                        }
                } else {
                        if( out(i, g) < (1-maf(g))*(1-maf(g))/4 ){
                                out(i, g) = 0;
                        } else if( out(i, g) < (1-maf(g))*(1-maf(g))*3/4 + maf(g)*(1-maf(g)) ){
                                out(i, g) = 1;
                        } else {
                                out(i, g) = 2;
                        }
                }
            }
        }
        return(out);
}", depends="RcppArmadillo")

calculateEffectsG = function(GenomeMatrix, BetaMatrix){
        ## takes SNP dosages for each individual and a matrix of betas for each SNP and TIME
        ## and multiplies them to get a matrix of INDxTIME
        EffectsG = GenomeMatrix %*% as.matrix(BetaMatrix)

        return(EffectsG)
}

## produces environmental effects for each day
## based on Gompertz distribution
cppFunction('arma::mat sumEffects(arma::mat genome, arma::vec env){
        arma::mat envT = env.t();
        for(int i=0; i<genome.n_rows; i++){
                genome.row(i) += envT;
        }
        return(genome);
}', depends="RcppArmadillo")

## a bit faster call to exp()
cppFunction('arma::vec exp_own(arma::vec min){
        return(exp(min));
}', depends="RcppArmadillo")

getSurvivals = function(EffectsT, lambda){
        ## takes the known effects and introduces stochasticity
        ## i.e. tests the e^hazards against Unif(0,1)
        risk = 1:nrow(EffectsT)
        time = 1
        SurvivalTimes = rep(ncol(EffectsT), length(risk))
        
        ## for each individual, test if it survives another day
        while(length(risk)>0 & time<=ncol(EffectsT)){
                p = exp_own(-lambda * exp_own(EffectsT[risk,time]))
                live = runif(length(p)) <= p
                SurvivalTimes[risk[!live]] = time
                risk = risk[live]
                time = time + 1
        }
        
        return(SurvivalTimes)
}

############ simulate unrelated individuals

simulateWithEffect = function(SnpTable){
        SnpMatrix = snpTableToMatrix(SnpTable)
        BetaMatrix = snpMatrixToBetas(SnpMatrix, maxTime)
        
        GenomeMatrix = generateGenome(inds, SnpMatrix$maf)
        
        EffectsE = alpha * 1:(maxTime-150)
        EffectsG = calculateEffectsG(GenomeMatrix, BetaMatrix)
        EffectsT = sumEffects(EffectsG, EffectsE)
        survTimes = getSurvivals(EffectsT, lambda)
        
        return(data.frame(GA=survTimes)+150)
}


############ simulate sibs

simulateWithEffectSibs = function(SnpTable){
        SnpMatrix = snpTableToMatrix(SnpTable)
        BetaMatrix = snpMatrixToBetas(SnpMatrix, maxTime)

        GenomeMatrix = generateGenome(inds, SnpMatrix$maf)
        GenomeMatrixSib = generateSibGenome(GenomeMatrix, SnpMatrix$maf)
        
        EffectsE = alpha * 1:(maxTime-150)
        EffectsG = calculateEffectsG(GenomeMatrix, BetaMatrix)
        EffectsT = sumEffects(EffectsG, EffectsE)
        survTimes = getSurvivals(EffectsT, lambda)

        EffectsG = calculateEffectsG(GenomeMatrixSib, BetaMatrix)
        EffectsT = sumEffects(EffectsG, EffectsE)
        survTimes2 = getSurvivals(EffectsT, lambda)
        return(data.frame(GA.x=survTimes, GA.y=survTimes2)+150)
}

calculateQuantiles = function(SurvTimes){
        ## calculates whatever quantiles are needed for the heritability estimation
        SurvSummary = mutate(SurvTimes, ga1_bin = GA.x %/% 7 * 7) %>%
                filter(ga1_bin>=230, ga1_bin<=290) %>%
                group_by(ga1_bin) %>%
                summarize(Q5 = quantile(GA.y, 0.05), Q25 = quantile(GA.y, 0.25), Q50 = quantile(GA.y, 0.50),
                          Q75 = quantile(GA.y, 0.75), Q95 = quantile(GA.y, 0.95), n=n())
        return(SurvSummary)
}

## user sets censoring time
maxTime = 310
## user sets sample size (will be multiplied by 12)
inds = 3e4
## alpha fitted from swed MFR
alpha = 0.12124
## lambda fitted from swed MFR
lambda = exp(-18.395) ## using 150 d cutoff, adjusted MFR. otherwise 6.05e-18

### set path to the quantiles calculated from observed data
RealSummary = read.table()

###############
## this launcher is to be used for automatic optimization
fitModelSibs = function(SnpTable){
        tt = mclapply(1:12, function(x) simulateWithEffectSibs(SnpTable), mc.cores=6)
        SurvTimes = data.frame(GA.x = c(sapply(tt, "[[", 1)), GA.y = c(sapply(tt, "[[", 2)))

        SurvSummary = calculateQuantiles(SurvTimes)
        ModelSummary = inner_join(SurvSummary, RealSummary, by="ga1_bin") %>%
                mutate(sq5 = (Q5.x-Q5.y)^2*sqrt(n.y), sq25 = (Q25.x-Q25.y)^2*sqrt(n.y), sq50 = (Q50.x-Q50.y)^2*sqrt(n.y),
                       sq75 = (Q75.x-Q75.y)^2*sqrt(n.y), sq95 = (Q95.x-Q95.y)^2*sqrt(n.y)) %>%
                mutate(ssq_bin = (sq5 + sq25 + sq50 + sq75 + sq95)/sqrt(n.y))
        costs = c(sum(ModelSummary$sq5), sum(ModelSummary$sq25), sum(ModelSummary$sq50),
                  sum(ModelSummary$sq75), sum(ModelSummary$sq95))/nrow(ModelSummary)
        cost = sum(costs)
        
        p1 = ggplot(ModelSummary, aes(x=ga1_bin)) + coord_cartesian(xlim=c(230, 290), ylim=c(220,300)) + theme_gray() +
                ggtitle(sprintf("Sum of Squares: %.5g \nEffect: %g / %g", cost, SnpTable$effect[2], SnpTable$effect[3])) +
                geom_line(aes(y=Q5.y)) + geom_line(aes(y=Q25.y)) + geom_line(aes(y=Q50.y)) +
                geom_line(aes(y=Q75.y)) + geom_line(aes(y=Q95.y)) +
                geom_line(aes(y=Q5.x), col="turquoise", size=1) + geom_line(aes(y=Q25.x), col="turquoise", size=1) + 
                geom_line(aes(y=Q50.x), col="turquoise", size=1) +  geom_line(aes(y=Q75.x), col="turquoise", size=1) + 
                geom_line(aes(y=Q95.x), col="turquoise", size=1) + 
                geom_bar(stat="identity", aes(y=230+50*ssq_bin/(1000+ssq_bin)), fill="coral", width=2)
        print(p1)
        
        return(cost)
}

## usage:
SnpTable = data.frame(name=c("filler", "fixed", "varying"), effect=c(0, 1, 0),
                      spread=c(1, 1, 10), number=c(1, 6, 1),
                      maf=c(0.5, 0.5, 0.5), time=c(50, 250-150, 250-150),
                      type=c("uniform", "uniform", "normal"))
fitModelSibs(SnpTable)

simulateWithEffect(SnpTable)
