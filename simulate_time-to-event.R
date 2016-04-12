library(dplyr)
library(ggplot2)
library(cowplot)
library(parallel)


simulator = function(timeframe, sampleSize, lambda, alpha, effect, timeOn, timeOff){
        def = data.frame(time = timeframe, hazard = 0, cumhazard = 0, survival = 1)
        n = 1
        for(t in def$time){
                ## generate hazard and cumulative hazard functions
                ## (at the moment both are generated analitically, but possible to just add hazards to get the cumulative hazard)
                if(t<timeOn){
                        def$hazard[n] = lambda * exp(alpha*t)
                        def$cumhazard[n] = lambda / alpha * (exp(alpha*t) - 1)
                } else if (t>=timeOn & t<timeOff){
                        def$hazard[n] = lambda * exp(alpha*t + effect)
                        def$cumhazard[n] = lambda / alpha * (exp(alpha*timeOn) - 1 + exp(effect + alpha*t) - exp(effect + alpha*timeOn))
                } else if (t>=timeOff){
                        def$hazard[n] = lambda * exp(alpha*t)
                        def$cumhazard[n] = lambda / alpha * (exp(alpha*timeOn) - 1 + exp(effect + alpha*timeOff) - 
                                                       exp(effect + alpha*timeOn) + exp(alpha*t) - exp(alpha*timeOff))
                }
                ## calculate survival function
                def$survival[n] = exp(-def$cumhazard[n])
                n = n + 1
        }
        ## plot the defined functions
        plot1 = ggplot(def, aes(x=time)) + geom_line(aes(y=hazard), col="red") + geom_line(aes(y=cumhazard), col="darkorange") + 
                geom_line(aes(y=survival), col="green") + 
                coord_cartesian(xlim=c(200,300))
        
        ## define range upper limits and inverse cum. hazards for each range
        range1 = def$cumhazard[timeOn]
        range2 = def$cumhazard[timeOff]
        inverse1 = function(t){
                return(log(t * alpha / lambda + 1) / alpha)
        }
        inverse2 = function(t){
                temp = t * alpha / (lambda * exp(effect)) - exp(alpha*timeOn - effect) + exp(alpha*timeOn) + exp(-effect)
                return(log(temp) / alpha)
        }
        inverse3 = function(t){
                temp = t * alpha / lambda - exp(alpha*timeOn) - exp(alpha*timeOff + effect) + exp(alpha*timeOn + effect) +
                        exp(alpha*timeOff) + 1
                return(log(temp) / alpha)
        }
        
        ## simulate event time distribution
        t = -log(runif(sampleSize))
        events = NULL
        for(i in 1:length(t)){
                if(t[i]<range1){
                        events[i] = inverse1(t[i])
                } else if (t[i]<range2){
                        events[i] = inverse2(t[i])
                } else {
                        events[i] = inverse3(t[i])
                }
        }
        plot2 = qplot(events, binwidth=1) + coord_cartesian(xlim=c(200,300))
        plot_grid(plot1, plot2, nrow=2)
}

simulator(1:300, 10000, 6.05e-18, 0.133, 0.5, 240, 285)

#####################################################
#####################################################
## binomial simulator for any empirical h(t)

snpTableToMatrix = function(SnpTable){
        ## convert a table with numbers of SNPs of different classes
        ## and output a matrix with the effects of individual SNPs
        SnpMatrix = data.frame(name = rep(SnpTable$name, SnpTable$number),
                            effect = rep(SnpTable$effect, SnpTable$number),
                            spread = rep(SnpTable$spread, SnpTable$number),
                            time = rep(SnpTable$time, SnpTable$number))

        return(SnpMatrix)
}

snpMatrixToBetas = function(SnpMatrix, maxTime, type){
        ## calculates beta coefficients for each position in time
        ## based on the spread on timing of each SNPs effect
        BetaMatrix = matrix(0, nrow=nrow(SnpMatrix), ncol=maxTime-150)
        if(type=="normal"){
                for(s in 1:nrow(SnpMatrix)){
                        m = SnpMatrix$time[s]
                        sd = SnpMatrix$spread[s]
                        b = SnpMatrix$effect[s]
                        BetaMatrix[s,] = dnorm(1:(maxTime-150), mean=m, sd=sd)*b
                }
        } else if(type=="uniform"){
                for(s in 1:nrow(SnpMatrix)){
                        b = SnpMatrix$effect[s]
                        BetaMatrix[s,] = rep(b, maxTime-150)
                }
        }
        BetaMatrix[BetaMatrix<1e-3] = 0
        BetaMatrix = Matrix(BetaMatrix)
        return(BetaMatrix)        
}

simulateGenome = function(inds, mafs){
        ## takes number of individuals, MAFs (SNP number estimated from MAF vector length)
        ## generates SNPs for each individual based only on MAFs, no LD
        GenomeMatrix = sapply(mafs, function(x) sample(c(0,1,2), inds, replace = TRUE, prob = c((1-x)^2, x*(1-x), x^2) ))
        GenomeMatrix = Matrix(GenomeMatrix)
        
        return(GenomeMatrix)
}

generateSibGenome = function(snpVector, maf){
        ## MAF is q (allele 0 frequency)
        zeros = snpVector==0
        ones = snpVector==1
        twos = snpVector==2
        snpVector[zeros] = sample(c(0,1), sum(zeros), replace=TRUE, prob=c(1-maf, maf))
        snpVector[ones] = sample(c(0,1,2), sum(ones), replace=TRUE, prob=c((1-maf)/2, 1/2, maf/2))
        snpVector[twos] = sample(c(1,2), sum(twos), replace=TRUE, prob=c(1-maf, maf))
        
        return(snpVector)
}

calculateEffectsG = function(GenomeMatrix, BetaMatrix){
        ## takes SNP dosages for each individual and a matrix of betas for each SNP and time
        ## and multiplies them to get a matrix of INDxTIME
        EffectsG = GenomeMatrix %*% BetaMatrix
        
        return(EffectsG)
}

calculateEffectsE = function(maxTime, alpha, inds){
        ## produces environmental effects for each day
        ## based on Gompertz distribution
        EffectsE = matrix(rep(alpha * 1:(maxTime-150), inds), nrow=inds, byrow=TRUE)

        return(EffectsE)
}

sumEffects = function(EffectsG, EffectsE){
        ## just add
        EffectsT = as.matrix(EffectsG) + EffectsE
        return(EffectsT)
}

getSurvivals = function(EffectsT, lambda){
        ## takes the known effects and introduces stochasticity
        ## i.e. tests the e^hazards against Unif(0,1)
        risk = 1:nrow(EffectsT)
        time = 1
        SurvivalTimes = rep(ncol(EffectsT), length(risk))
        
        ## for each individual, test if it survives another day
        while(length(risk)>0 & time<=ncol(EffectsT)){
                p = exp(-lambda * exp(EffectsT[risk,time]))
                live = runif(length(p)) <= p
                SurvivalTimes[risk[!live]] = time
                risk = risk[live]
                time = time + 1
        }

        return(SurvivalTimes)
}

## user sets censoring time
maxTime = 310
## user sets sample size
inds = 1e4
## alpha fitted from swed MFR
alpha = 0.133
## lambda fitted from swed MFR
lambda = 2.69e-09 ## using 150 d cutoff. otherwise 6.05e-18

simulateWithEffect = function(arg_eff, arg_time){
        ## user sets this table
        SnpTable = data.frame(name=c("early", "medium", "late"), effect=c(0, 0, arg_eff),
                              spread=c(1,1,1), number=c(1,1,1), time=c(150, 200, arg_time))
        SnpMatrix = snpTableToMatrix(SnpTable)
        BetaMatrix = snpMatrixToBetas(SnpMatrix, maxTime, "normal")
        ## use runif for MAF later
        GenomeMatrix = simulateGenome(inds, rep(1, nrow(SnpMatrix)))
        EffectsG = calculateEffectsG(GenomeMatrix, BetaMatrix)
        EffectsE = calculateEffectsE(maxTime, alpha, inds)
        EffectsT = sumEffects(EffectsG, EffectsE)
        survTimes = getSurvivals(EffectsT, lambda)
        return(survTimes)
#         return(data.frame(eff=arg_eff,time=arg_time, m=mean(survTimes), s=sd(survTimes),
#                           q10=quantile(survTimes, 0.1), q90=quantile(survTimes, 0.9)))
}

eff = c(0, 5, 10, 14, 17, 20, 22, 24, 26, 28, 30)
times = c(160, 180, 200, 220, 240, 260, 280)
arg_df = NULL
for(e in eff){
        for (t in times){
                arg_df = rbind(arg_df, data.frame(e, t))
        }
}

rez = mcmapply(simulateWithEffect, arg_df$e, arg_df$t, mc.cores = 6)
rez = t(rez)
rez = data.frame(rez)
rez2 = data.frame(eff = unlist(rez$eff),
                  time = unlist(rez$time),
                  m = unlist(rez$m), s = unlist(rez$s),
                  q10 = unlist(rez$q10), q90 = unlist(rez$q90))

ggplot(rez2, aes(x=eff)) + geom_line(aes(y=m)) + geom_ribbon(aes(ymin=m-1.96*sd, ymax=m+1.96*sd), fill="grey", alpha=0.7) + facet_grid(.~time) + theme_grey()
ggplot(rez2, aes(x=eff, group=time)) + geom_line(aes(y=m)) + geom_ribbon(aes(ymin=m-1.96*sd, ymax=m+1.96*sd, fill=factor(time)), alpha=0.6) +
        scale_fill_brewer() + theme_grey()
ggplot(rez2[rez2$eff>0,], aes(x=time, col=eff)) + geom_point(aes(y=(m-278.6)/exp(eff)))
mutate(rez2, beta=(m-278.6)/exp(eff))

qplot(survTimes, binwidth=1)
mean(survTimes); sd(survTimes)

paramSimulate = function(n, alpha, lambda){
        return(data.frame(dat = log(-log(runif(n)) * alpha / lambda + 1) / alpha,
                          alpha = alpha, lambda = lambda))
}

simmean = 0.1326341
simrate = exp(-19.73372)
simdat = bind_rows(est = paramSimulate(nrow(mfr), simmean, simrate),
                   low = paramSimulate(nrow(mfr), simmean-1.96*9.28e-05, simrate-1.96*3.23e-11),
                   upp = paramSimulate(nrow(mfr), simmean+1.96*9.28e-05, simrate+1.96*3.23e-11), .id="run")
ggplot(data=simdat, aes(x=dat+150, group=factor(run))) +geom_line(stat="bin", binwidth=1, aes(col=factor(run))) + scale_x_continuous(limits=c(220,250))
group_by(simdat, run) %>% summarize(mean(dat), sd(dat))

##

eff0 = simulateWithEffect(0, 100)+150
eff1 = simulateWithEffect(1, 100)+150
eff1 = simulateWithEffect(2, 100)+150
ggplot() + geom_histogram(aes(x=eff0), col="red", alpha=0.1, binwidth=1) +
        geom_histogram(aes(x=eff1), col="blue", alpha=0.1, binwidth=1) +
        geom_histogram(aes(x=eff1), col="pink", alpha=0.1, binwidth=1)

########## optimizer part

fw <- function(x){
        sim = simulateWithEffect(x, 200)
        abs(sim$m-240)
}

res <- optim(1, fw, method = "BFGS", control=list(trace=1, REPORT=1))
res
?optim

############ simulate sibs

simulateWithEffectSibs = function(arg_eff, arg_time, arg_num, arg_spread, arg_maf, arg_mode){
        ## user sets this table
        SnpTable = data.frame(name=c("early", "medium", "late"), effect=c(0, 0, arg_eff),
                              spread=c(1,1,arg_spread), number=c(1,1,arg_num), time=c(150, 200, arg_time))
        SnpMatrix = snpTableToMatrix(SnpTable)          ## 1.1 ms
        BetaMatrix = snpMatrixToBetas(SnpMatrix, maxTime, arg_mode)     ## 0.4 ms
        ## use runif for MAF later
        GenomeMatrix = simulateGenome(inds, rep(arg_maf, nrow(SnpMatrix)))      ## 3.5 ms
        GenomeMatrixSib = apply(as.matrix(GenomeMatrix), 2, generateSibGenome, arg_maf)         ## 4 ms
        
        EffectsG = calculateEffectsG(GenomeMatrix, BetaMatrix)          ## 10 ms
        EffectsE = calculateEffectsE(maxTime, alpha, inds)              ## 19 ms
        EffectsT = sumEffects(EffectsG, EffectsE)             ## 60 ms
        survTimes = getSurvivals(EffectsT, lambda)      ## 76 ms
        
        EffectsG = calculateEffectsG(GenomeMatrixSib, BetaMatrix)       ## 8 ms
        EffectsE = calculateEffectsE(maxTime, alpha, inds)      ## 18 ms
        EffectsT = sumEffects(EffectsG, EffectsE)       ## 60 ms
        survTimes2 = getSurvivals(EffectsT, lambda)
        return(data.frame(GA.x=survTimes, GA.y=survTimes2))
}

calculateQuantiles = function(SurvTimes){
        ## calculates whatever quantiles are needed for the heritability estimation
        SurvSummary = mutate(SurvTimes, ga1_bin = GA.x %/% 10 * 10) %>%
                group_by(ga1_bin) %>%
                summarize(Q5 = quantile(GA.y, 0.05), Q50 = quantile(GA.y, 0.50), Q95 = quantile(GA.y, 0.95))
        return(SurvSummary)
}

calcCost = function(ModelSummary, RealSummary){
        ## cost function: sum-of-squares of three quantiles
        
        
        return(cost)
}

RealSummary = calculateQuantiles(cousins_mean4)

palette = colorRampPalette(colors=c("green", "yellow", "red"))(10)

## initial settings and placeholders
inds = 5e4
bestcosts = c(1e5, 1e5, 1e5)
bestcost = 3e5
allCosts = NULL
params = seq(0, 3, by=0.1)

for(i in params){
        SurvTimes = simulateWithEffectSibs(i, 100, 1, 1, 0.5, "uniform")+150
        qplot(SurvTimes$GA.x)
        SurvSummary = calculateQuantiles(SurvTimes)
        ModelSummary = inner_join(SurvSummary, RealSummary, by="ga1_bin") %>%
                mutate(sq5 = (Q5.x-Q5.y)^2, sq50 = (Q50.x-Q50.y)^2, sq95 = (Q95.x-Q95.y)^2)
        costs = c(sum(ModelSummary$sq5), sum(ModelSummary$sq50), sum(ModelSummary$sq95))/nrow(ModelSummary)
        cost = sum(costs)
        allCosts = c(allCosts, cost)
        
        if(cost<bestcost){
                ## assign parameters as best
                bestcosts = costs
                bestcost = cost
        }

        ## plot, for fun
        colors = palette[ round(costs/bestcosts/10)+1 ]
        p1 = ggplot(ModelSummary, aes(x=ga1_bin)) + coord_cartesian(xlim=c(200, 290)) + theme_gray() +
                ggtitle(sprintf("Sum of Squares: %.5g \nEffect: %g", cost, i)) +
                geom_line(aes(y=Q5.y)) + geom_line(aes(y=Q50.y)) + geom_line(aes(y=Q95.y)) +
                geom_line(aes(y=Q5.x), col=colors[1]) + geom_line(aes(y=Q50.x), col=colors[2]) + geom_line(aes(y=Q95.x), col=colors[3])
        print(p1)
}

## plot goodness-of-fit
ggplot(data.frame(params, allCosts), aes(x=params)) + geom_line(aes(y=allCosts))



## plot inheritance curve
regrline = lm(GA.y ~ GA.x, data=sibs_point)
ggplot(rez, aes(x=ga1_bin)) + geom_line(aes(y=MeanGA), size=2) + geom_line(aes(y=Q5)) + geom_line(aes(y=Q95)) +
        coord_cartesian(xlim=c(200, 290)) +
        geom_abline(slope=regrline$coefficients[2], intercept=regrline$coefficients[1], col="red")


temp = cousins_mean2
dim(temp)
temp$ga1_bin = temp$GA.x %/% 10 * 10
sort(unique(temp$ga1_bin))

rez = group_by(temp, ga1_bin) %>%  summarize(N = n(), MeanGA = mean(GA.y),
                                                   Q5 = quantile(GA.y, 0.05), Q50 = quantile(GA.y, 0.50), Q95 = quantile(GA.y, 0.95))
rez
regrline = lm(GA.y ~ GA.x, data=temp)
ggplot(rez, aes(x=ga1_bin)) + geom_line(aes(y=MeanGA), size=2) + geom_line(aes(y=Q5)) + geom_line(aes(y=Q95)) +
        coord_cartesian(xlim=c(200, 290)) +
        geom_abline(slope=regrline$coefficients[2], intercept=regrline$coefficients[1], col="red")
