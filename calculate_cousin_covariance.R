options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(cowplot)
library(parallel)
library(mvtnorm)

### this script reads in the swed mfr file, pairs of full-sib cousins,
### and calculates covariance across different gestational ages

cousins = read.table("/home/lab/TRANSPORTER_JULIUI/cousins_pairs_maternal.csv", h=T)

head(cousins); dim(cousins) # 240 k pairs

## here we see that some half sibs are left
c1 = select(cousins, matches("fid|1$"))
c2 = select(cousins, matches("fid|2$"))
names(c2) = names(c1)
c1 = bind_rows(c1, c2) # 480 k total
group_by(c1, mom1) %>%
        summarize(ndads = length(unique(dad1))) %>%
        arrange(desc(ndads))

## remove half sibs, leaving only the children from the first father
cousins_firstdad = arrange(c1, mom1, birthyear1) %>%
        group_by(mom1) %>%
        mutate(dad_rank = rank(dad1, ties.method="min")) %>%
        filter(dad_rank==1)    ## 446 k remain
## remove duplicate kid ids
cousins_firstdad = cousins_firstdad[which(!duplicated(cousins_firstdad$kid1)),] ## 198 k remain (211 k uniques in the beginning)


### mfr was extracted from Verena's mfr using this query:
#       select p.lpnr_barn, LopNrMor, LopNrFar, GRDBS, p.KON, PARITET, MLANGD, MVIKT, MALDER, FAMSIT, MFODLAND, ROK0, ROK1, ROK2
#       from mfr15 m INNER JOIN parents p ON p.lpnr_barn=m.lpnr_BARN 
#       WHERE FLSPONT=1 AND MISSB=0 AND GRDBS>0 AND BORDF2=1 AND DODFOD=0 INTO OUTFILE '/tmp/output.csv';

mfr = read.table("~/Documents/swed_mfr_games/mfr_verena_cousins_spont.csv", h=T, fill=T, sep="\t")
head(mfr); dim(mfr)

summarizeWindow = function(df, window){
        out = NULL
        df$GA_m = round((df$GA.x + df$GA.y)/2)
        df = select(df, starts_with("GA"))
        
        for(w in -window:window){
                df$GA_w = w
                out = bind_rows(out, mutate(df, GA_gr=GA_m + GA_w))
        }
        out = filter(out, GA_gr>=min(GA_m), GA_gr<=max(GA_m)) %>%
                group_by(GA_gr) %>%
                summarize(c = cov(GA.x, GA.y), n = n()) %>%
                mutate(ci_low = c*(n-1)/qchisq(0.975, n-1), ci_upp = c*(n-1)/qchisq(0.025, n-1)) %>%
                filter(!is.na(c))
        return(out)
}

## attach GA to the cousins
cousins = inner_join(cousins_firstdad, mfr, by=c("kid1" = "id_barn")) ## 151 k remain, most losses due to FLSPONT=0 or FLINDUKT=1
cousins_mean = group_by(cousins, mom1) %>%
        summarize(fid = max(fid), GA = round(mean(GA)))   ## 95 k unique mothers remain
## form sister pairs
cousins_mean = inner_join(cousins_mean, cousins_mean, by="fid") %>%
        filter(mom1.x>mom1.y)     ## 52 k families, but not all have two sisters, so 47 k remain
## if there is more than one sister, randomly grab one
cousins_mean = cousins_mean[!duplicated(cousins_mean$fid),]    ## 39 k remain
cov(cousins_mean$GA.x, cousins_mean$GA.y) ## = 1/8 V_F + 1/2 V_M; V_M = 0.2297
cov(cousins_mean$GA.x, cousins_mean$GA.y)/var(c(cousins_mean$GA.x, cousins_mean$GA.y))*8/5
        
        
#### option b) taking the child with lowest parity
cousins_mean2 = group_by(cousins, mom1) %>%
        top_n(1, desc(PARITET)) %>% ungroup()
cousins_mean2 = inner_join(cousins_mean2, cousins_mean2, by="fid") %>%
        filter(mom1.x>mom1.y)
cousins_mean2 = cousins_mean2[!duplicated(cousins_mean2$fid),]
cov(cousins_mean2$GA.x, cousins_mean2$GA.y) ## = 1/8 V_F + 1/2 V_M; V_M = 0.2244
cov(cousins_mean2$GA.x, cousins_mean2$GA.y)/var(c(cousins_mean2$GA.x, cousins_mean2$GA.y))*8/5


#### option c) adjusting for covariates + mean within mother
mfrc = filter(mfr, PARITET>0, ROK1>0, MHEIGHT>143, MWEIGHT>40 & MWEIGHT<150)
mfrc$MFODLAND_cat = mfrc$MFODLAND!="SVERIGE"
mfrc$FAMSIT_cat = mfrc$FAMSIT!=1
mfrc$MAGE_cat = cut(mfrc$MAGE, breaks=c(0,20, 25, 30, 35, 100), labels = 1:5)
mfrc$MAGE_cat = factor(mfrc$MAGE_cat, levels = c(2,1,3,4,5))
mfrc$PARITET[mfrc$PARITET>4] = 4

lm_formula = as.formula("GA ~ factor(PARITET) + KON + MHEIGHT + YEAR + MWEIGHT + FAMSIT_cat + MAGE_cat + MFODLAND_cat + factor(ROK1)")

summary(lm(lm_formula, data=mfrc))

mfr_fit = lm(lm_formula, data=mfrc)$residuals + mean(mfrc$GA)
mfr_fit = data.frame(id_barn = as.integer(mfrc$id_barn), GA = mfr_fit, PARITET = mfrc$PARITET)

cousins = inner_join(cousins_firstdad, mfr_fit, by=c("kid1" = "id_barn")) ## 136 k remain
cousins_mean3 = group_by(cousins, mom1) %>%
        summarize(fid = max(fid), GA = round(mean(GA)))      ## 89 k reamin
cousins_mean3 = inner_join(cousins_mean3, cousins_mean3, by="fid") %>%
        filter(mom1.x>mom1.y)      ## 42 k remain
cousins_mean3 = cousins_mean3[!duplicated(cousins_mean3$fid),]  ## 36 k remain
cov(cousins_mean3$GA.x, cousins_mean3$GA.y) ## = 1/8 V_F + 1/2 V_M; V_M = 0.2207
cov(cousins_mean3$GA.x, cousins_mean3$GA.y)/var(c(cousins_mean3$GA.x, cousins_mean3$GA.y))*8/5


### option d) adjusting for covariates + first born

cousins_mean4 = group_by(cousins, mom1) %>%
        top_n(1, desc(PARITET)) %>% ungroup()
cousins_mean4 = inner_join(cousins_mean4, cousins_mean4, by="fid") %>%
        filter(mom1.x>mom1.y)
cousins_mean4 = cousins_mean4[!duplicated(cousins_mean4$fid),]
cov(cousins_mean4$GA.x, cousins_mean4$GA.y) ## = 1/8 V_F + 1/2 V_M; V_M = 0.2160
cov(cousins_mean4$GA.x, cousins_mean4$GA.y)/var(c(cousins_mean4$GA.x, cousins_mean4$GA.y))*8/5

cousins_mean_sum = calculateQuantiles(cousins_mean)
cousins_mean_sum2 = calculateQuantiles(cousins_mean2)
cousins_mean_sum3 = calculateQuantiles(cousins_mean3)
cousins_mean_sum4 = calculateQuantiles(cousins_mean4)
write.table(cousins_mean_sum, "~/Documents/swed_mfr_games/simulations_local/finalrez/cousins_mean_sum.csv", quote = F, row.names = F, col.names = T)
write.table(cousins_mean_sum2, "~/Documents/swed_mfr_games/simulations_local/finalrez/cousins_mean_sum2.csv", quote = F, row.names = F, col.names = T)
write.table(cousins_mean_sum3, "~/Documents/swed_mfr_games/simulations_local/finalrez/cousins_mean_sum3.csv", quote = F, row.names = F, col.names = T)
write.table(cousins_mean_sum4, "~/Documents/swed_mfr_games/simulations_local/finalrez/cousins_mean_sum4.csv", quote = F, row.names = F, col.names = T)
write.table(tidy(lmfit), "~/Documents/swed_mfr_games/simulations_local/finalrez/mfr_fit.csv", quote = F, row.names = F, col.names = T)

p1 = ggplot(cousins_mean_sum, aes(x=GA_gr)) + theme(panel.grid.major = element_line(linetype="dashed", color="grey50", size=0.5)) +
        geom_ribbon(aes(ymax=-ci_upp, ymin=-ci_low), fill="grey70", alpha=0.4) + coord_cartesian(ylim=c(-50, 4000)) +
        geom_line(aes(y=-c), col="black") + ylab("-covariance")
p2 = ggplot(cousins_mean_sum2, aes(x=GA_gr)) + theme(panel.grid.major = element_line(linetype="dashed", color="grey50", size=0.5)) +
        geom_ribbon(aes(ymax=-ci_upp, ymin=-ci_low), fill="grey70", alpha=0.4) + coord_cartesian(ylim=c(-50, 4000)) +
        geom_line(aes(y=-c), col="black") + ylab("-covariance")
p3 = ggplot(cousins_mean_sum3, aes(x=GA_gr)) + theme(panel.grid.major = element_line(linetype="dashed", color="grey50", size=0.5)) +
        geom_ribbon(aes(ymax=-ci_upp, ymin=-ci_low), fill="grey70", alpha=0.4) + coord_cartesian(ylim=c(-50, 4000)) +
        geom_line(aes(y=-c), col="black") + ylab("-covariance")
p4 = ggplot(cousins_mean_sum4, aes(x=GA_gr)) + theme(panel.grid.major = element_line(linetype="dashed", color="grey50", size=0.5)) +
        geom_ribbon(aes(ymax=-ci_upp, ymin=-ci_low), fill="grey70", alpha=0.4) + coord_cartesian(ylim=c(-50, 4000)) +
        geom_line(aes(y=-c), col="black") + ylab("-covariance")

#### all plots can be nicely seen here
plot_grid(p1+ggtitle("unadjusted, mean of all children")+geom_hline(yintercept=0, col="red", size=0.2),
          p2+ggtitle("unadjusted, only first child")+geom_hline(yintercept=0, col="red", size=0.2),
          p3+ggtitle("adjusted, mean of all children")+geom_hline(yintercept=0, col="red", size=0.2),
          p4+ggtitle("adjusted, only first child")+geom_hline(yintercept=0, col="red", size=0.2))

#### now do everything with shuffling
shuffleCousinTimes = function(df, dfsum, window, nperm){
        pb = txtProgressBar(1, nperm)
        
        for(i in 1:nperm){
                setTxtProgressBar(pb, i)
                df$GA.x = sample(df$GA.x, nrow(df))
                df$GA.y = sample(df$GA.y, nrow(df))
                df$GA_m = round((df$GA.x+df$GA.y)/2)
                
                loop_sum = summarizeWindow(df, window) %>% select(one_of("GA_gr","c"))
                names(loop_sum)[2] = paste("cs", i, sep="")
                dfsum = left_join(dfsum, loop_sum, by="GA_gr")
        }
        
        dfsum = dfsum[which(apply(dfsum[,6:(nperm+5)], 1, function(x) !all(is.na(x)))), ]
        dfsum$shuffquantL = apply(dfsum[,6:(nperm+5)], 1, quantile, 0.025, na.rm=T)
        dfsum$shuffquantU = apply(dfsum[,6:(nperm+5)], 1, quantile, 0.975, na.rm=T)
        dfsum$shuffmean = apply(dfsum[,6:(nperm+5)], 1, mean, na.rm=T)
        dfsum = mutate(dfsum, c_ratio = shuffmean/c, c_dif = c-shuffmean)
        dfsum = select(dfsum, -starts_with("cs"))
        
        return(dfsum)
}

shuffResults = mclapply(list(m1=list(df=cousins_mean, dfs=cousins_mean_sum),
            m2=list(df=cousins_mean2, dfs=cousins_mean_sum2),
            m3=list(df=cousins_mean3, dfs=cousins_mean_sum3),
            m4=list(df=cousins_mean4, dfs=cousins_mean_sum4)), function(x) shuffleCousinTimes(x$df, x$dfs, 7, 100), mc.cores=4)

shuffTimes = bind_rows(unadj_mean = shuffResults[[1]], unadj_first = shuffResults[[2]],
                       adj_mean = shuffResults[[3]], adj_first = shuffResults[[4]], .id="model")

### all plots showing shuffled means and CIs
ggplot(shuffTimes, aes(x=GA_gr)) + theme(panel.grid.major = element_line(linetype="dashed", color="grey50", size=0.5)) +
        geom_ribbon(aes(ymax=-shuffquantL, ymin=-shuffquantU), fill="skyblue", alpha=0.4) + 
        geom_line(aes(y=-c), col="black") + geom_line(aes(y=-shuffmean), col="blue") +
        facet_wrap(~model) + 
        coord_cartesian(ylim=c(-10, 4000)) + ylab("-covariance")

### comparing the ratios of all models
ggplot(shuffTimes, aes(x=GA_gr)) + geom_line(aes(y=c_ratio, col=model, alpha=as.numeric(shuffquantU<c | shuffquantL>c))) +
        scale_alpha_continuous(guide="none", range = c(0.3, 1)) +
        coord_cartesian(ylim=c(0.1,2.5), xlim = c(200,300)) 

## comparing the control covariances of all models
ggplot(shuffTimes, aes(x=GA_gr)) + geom_line(aes(y=log(-shuffmean), col=model))


############ another way of introducing uniform covariance
### draw correlated mvnorm
### convert to quantiles
### take corresponding quantiles from cousin data
simCorrMath = function(dfin, window, nperm){
        sigma = cor(dfin[,c("GA.x", "GA.y")])
        sigma[c(2,3)] = sigma[c(2,3)]*1.15 ## bias correction
        #allGA = c(dfin$GA.x, dfin$GA.y)
        simnorm = rmvnorm(1e6, sigma=sigma)
        simnorm = pnorm(simnorm)
        sim = data.frame(GA.x = quantile(dfin$GA.x, simnorm[,1]), GA.y = quantile(dfin$GA.x, simnorm[,2]), row.names = NULL)
        #simsum = summarizeWindow(sim, window)
        #print(sigma); print(cov(sim)); print(cor(sim))
        #return(shuffleCousinTimes(sim, simsum, window, nperm))
        sim$GA.x = round(sim$GA.x); sim$GA.y = round(sim$GA.y)
        sim$GA.shuff = sample(sim$GA.y, nrow(sim))
        sim = group_by(sim, GA.x) %>% summarize(conc = sum(GA.y==GA.x), concshuff=sum(GA.shuff==GA.x), dens = n())
        return(sim)
}

simcorrs = mclapply(list(unadj_mean=cousins_mean, unadj_first=cousins_mean2,
                         adj_mean=cousins_mean3, adj_first=cousins_mean4), function(x) simCorrMath(x, 7, 100), mc.cores=4)
simcorrs2 = NULL
for(i in 1:length(simcorrs)){
        temp = simcorrs[[i]]
        temp$model = names(simcorrs)[i]
        simcorrs2 = bind_rows(simcorrs2, temp)
}

## confirmed: shuffled means produced from actual cousins and manually correlated cousin draws are equal everywhere >220
## meaning, we can compare the true means without caring about the shuffled ones
both = inner_join(shuffTimes, simcorrs2, by=c("model", "GA_gr"))
ggplot(both, aes(x=GA_gr, group=model)) + geom_line(aes(y=shuffmean.x/shuffmean.y, col=model)) + geom_hline(aes(yintercept=1), col="red")

ggplot(filter(both, model=="adj_first"), aes(x=GA_gr, group=model)) +
        geom_line(aes(y=c.y/c.x, col=model)) + geom_hline(aes(yintercept=1), col="red") + coord_cartesian(ylim=c(0.8,1.5))
ggplot(both, aes(x=GA_gr, group=model)) + geom_line(aes(y=c_ratio.x-c_ratio.y, col=model)) + geom_hline(aes(yintercept=0), col="red") +
        coord_cartesian(ylim=c(-1,1.5))

fit = lm(GA.y ~ GA.x, data=cousins_mean)
fit$ga_x = cousins_mean$GA.x

table(cousins_mean$GA.x == cousins_mean$GA.y) 


shuff = simCorrMath(cousins_mean, 0, 1)
shuff$concref = shuff$conc[shuff$GA.x==280]
shuff$concshuffref = shuff$concshuff[shuff$GA.x==280]
ggplot(shuff, aes(x=GA.x)) + geom_line(aes(y=conc/concshuff/concref*concshuffref), col="red")# + coord_cartesian(ylim=c(0, 0.0002))



#############################
###### transforming GAs into U(0,1)
summarizeWindowFixed = function(df, window, width){
        out = NULL
        df$GA_m = ((df$GA.x + df$GA.y)/2)
        df = select(df, starts_with("GA"))
        
        for(w in -window:window){
                df$GA_w = w
                out = bind_rows(out, mutate(df, GA_gr=round(GA_m + GA_w)))
        }
        out = filter(out, GA_gr>=min(GA_m)+window+width, GA_gr<=max(GA_m)-window-width,
                     abs(GA.x-GA_m)<=width, abs(GA.y-GA_m)<=width) %>%
                group_by(GA_gr) %>%
                summarize(c = sd(GA.x-GA.y), n = n()) %>%
                filter(!is.na(c))
        return(out)
}


out = NULL
for(i in 1:10){
        
        ## transform GAs
        unif1 = data.frame(GA.x = rank(cousins_mean2$GA.x, ties.method = "random"), GA.y = rank(cousins_mean2$GA.y, ties.method = "random"))
        unif1 = unif1 / nrow(unif1)
        cor(unif1)
        
        ## calculate cov of shuffled
        unif2 = mutate(unif1, GA.y=sample(GA.y, nrow(unif1))) * 100
        cor(unif2)
        unif2_sum = summarizeWindowFixed(unif2, 2, 0.5)
        # ggplot(unif2_sum, aes(x=GA_gr)) + geom_line(aes(y=-c))
        
        ## calculate cov with induced uniform correlation
        sigma = cor(unif1)
        sigma[2:3] = sigma[2:3]
        simnorm = rmvnorm(nrow(unif1), sigma=sigma)
        simnorm = pnorm(simnorm)
        #unif3 = data.frame(GA.x = quantile(unif1$GA.x, simnorm[,1]), GA.y = quantile(unif1$GA.x, simnorm[,2]), row.names = NULL)
        unif3 = data.frame(GA.x = simnorm[,1], GA.y = simnorm[,2])
        unif3 = unif3 * 100
        cor(unif3)
        summary(unif3)
        unif3_sum = summarizeWindowFixed(unif3, 2, 0.5)
        # ggplot(unif3_sum, aes(x=GA_gr)) + geom_line(aes(y=-c))
        
        ## calculate cov with induced point correlation
        unif4 = mutate(unif1, GA.y=sample(GA.y, nrow(unif1))) * 100
        unif4$GA.y = ifelse(rbinom(nrow(unif1), 1, 0.15)>0, unif4$GA.x, unif4$GA.y)
        cor(unif4)
        unif4_sum = summarizeWindowFixed(unif4, 2, 0.5)
        # ggplot(unif4_sum, aes(x=GA_gr)) + geom_line(aes(y=-c))
        
        ## calculate cov of actual data
        cor(unif1)
        unif1_sum = summarizeWindowFixed(unif1*100, 2, 0.5)
        # ggplot(unif1_sum, aes(x=GA_gr)) + geom_line(aes(y=-c))
        
        all = bind_rows(true=unif1_sum, shuffled = unif2_sum, corr_unif = unif3_sum, corr_unif_my = unif4_sum, .id="model")
        all = all[,1:3]
        if(i==1){
                out=all
        } else {
                out = left_join(out, all, by=c("model", "GA_gr"))
        }
}
out$cmean = apply(out[,3:12], 1, mean, na.rm=T)
out$q10 = apply(out[,3:12], 1, quantile, 0.025, na.rm=T)
out$q90 = apply(out[,3:12], 1, quantile, 0.975, na.rm=T)
out$trueGA = quantile(c(cousins_mean2$GA.x, cousins_mean2$GA.y), out$GA_gr/100)

ggplot(out, aes(x=GA_gr, group=model)) + geom_line(aes(y=cmean, color=model)) +
        geom_ribbon(aes(ymax=q10, ymin=q90, fill=model), alpha=0.2)

library(flexsurv)
mfr_gomp = mfr_fit
mfr_gomp_nocous = anti_join(mfr_gomp, cousins_mean4, by=c("id_barn"="kid1.x")) %>% anti_join(cousins_mean4, by=c("id_barn"="kid1.y"))
mfr_gomp_nocous$GA = mfr_gomp_nocous$GA-150
mfr_gomp_nocous$event = TRUE
gompfit = flexsurvreg(Surv(GA, event) ~ 1, data = mfr_gomp_nocous, dist = "gompertz")
gompfit$coefficients
gompfit$AIC

gompfit = flexsurvreg(Surv(GA, event) ~ 1, data = mfr_gomp_nocous, dist = "weibull")
gompfit$coefficients
gompfit$AIC

gompfit = flexsurvreg(Surv(GA, event) ~ 1, data = mfr_gomp_nocous, dist = "exp")
gompfit$coefficients
gompfit$AIC
