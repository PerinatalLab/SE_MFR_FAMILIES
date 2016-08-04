
library(tidyr)
library(corrplot)

SurvTimes = NULL
heights = seq(-2, 1, by=0.5)

for(h in heights){
        SnpTable = data.frame(name=c("filler", "fixed", "varying"), effect=c(0, 0, -h*15),
                              spread=c(1, 1, 10), number=c(1, 1, 1),
                              maf=c(0.5, 1, 1), time=c(50, 250-150, 250-150),
                              type=c("uniform", "uniform", "normal"))
        tt = mclapply(1:12, function(x) simulateWithEffect(SnpTable), mc.cores=6)
        tt = data.frame(GA = unlist(tt))
        tt$height = h
        SurvTimes = bind_rows(SurvTimes, tt)
}

ggplot(SurvTimes) + geom_density(aes(x=GA, col=height, group=height))

GAtable = mutate(SurvTimes, GAbin = cut(GA, 7)) %>% group_by(height, GAbin) %>% summarize(n=n()) %>% spread(GAbin, n)
GAtable = as.matrix(GAtable)
GAtable[is.na(GAtable)] = 0 
rownames(GAtable) = GAtable[,1]
GAtable = GAtable[,-1]
GAexp = as.numeric(table(SurvTimes$height)) %*% t(as.numeric(table(cut(SurvTimes$GA, 7)))) / nrow(SurvTimes)
GAenr = GAtable/GAexp
GAenr[GAtable<100] = 0

corrplot(GAenr, is.corr=F)
