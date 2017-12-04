library(ggplot2)

boxplot.thr = 0
hist.log.ymin = 10

hist = read.table('out/hotspot_histograms_100kb.txt')
names(hist) = c('pctl','sample','contacts','bin.count')
pctl.name = list('0'="Random",'90'="Top 10%", '95'="Top 5%", '99'="Top 1%")
hist$pctl = as.factor(hist$pctl)
hist$sample = paste(hist$pctl,hist$sample,sep='.')
hist$sample = as.factor(hist$sample)

bounds = data.frame(pctl=c(),x=c(),ymin=c(),ymax=c(),med=c())
box.data = data.frame(pctl=c(),contacts=c())

for (p in levels(hist$pctl)) {
    hp = hist[hist$pctl == p,]
    # Fill boundary/median data for plot
    for (x in unique(hp$contacts)) {
            bounds = rbind(bounds,data.frame(pctl=pctl.name[[toString(p)]],x=x,ymin=min(hp[hp$contacts == x,]$bin.count),ymax=max(hp[hp$contacts == x,]$bin.count),med=median(hp[hp$contacts == x,]$bin.count)))
    }
    # Fill boxplot data
    hp.thr = hp[hp$contacts >= boxplot.thr,]
    box.data = rbind(box.data, data.frame(pctl=pctl.name[[toString(p)]], contacts=rep(hp.thr$contacts,hp.thr$bin.count)))
}

# Plot all lines.
#ggplot(hist,aes(x=contacts,y=bin.count,color=pctl,group=sample)) + geom_line() + scale_y_log10() + xlim(c(0,50))

# Reorder factors
bounds$pctl = factor(bounds$pctl,c("Top 1%","Top 5%","Top 10%","Random"))
box.data$pctl = factor(box.data$pctl,c("Top 1%","Top 5%","Top 10%","Random"))

# Plot ribbon + median
pdf('hotspot_histograms.pdf',useDingbats=F,width=10,height=7)
ggplot(bounds) + geom_line(aes(x=x,y=med,color=pctl)) + geom_ribbon(aes(ymin=ymin,ymax=ymax,x=x,fill=pctl), alpha = 0.3) + xlim(c(0,20)) + scale_fill_discrete(name="Hotspots") + scale_color_discrete(name="Hotspots") + theme_minimal() + xlab('Hi-C contact frequency') + ylab('Pairwise interaction count')
dev.off()

# Boxplot for presentations (thr > 3)
pdf('hotspot_boxplot.pdf',useDingbats=F,width=7,height=10)
ggplot(box.data, aes(x=pctl,y=contacts,group=pctl,fill=pctl,color=pctl)) + geom_boxplot(outlier.shape=NA,alpha=0.3) + ylim(c(0,15)) + ylab('Hi-C contact frequency') +xlab("") + scale_fill_discrete(name="Hotspots") + scale_color_discrete(name="Hotspots") + theme_minimal() + theme(axis.text.x=element_blank())
dev.off()

# Log scale
bounds[bounds$ymax >= hist.log.ymin & bounds$ymin < hist.log.ymin,]$ymin = hist.log.ymin
pdf('hotspot_histograms_log.pdf',useDingbats=F,width=10,height=7)
ggplot(bounds) + geom_line(aes(x=x,y=med,color=pctl)) + geom_ribbon(aes(ymin=ymin,ymax=ymax,x=x,fill=pctl), alpha = 0.3) + xlim(c(0,25)) + scale_y_log10(limits=c(hist.log.ymin,max(bounds$ymax))) + scale_fill_discrete(name="Hotspots") + scale_color_discrete(name="Hotspots") + theme_minimal() + xlab('Hi-C contact frequency') + ylab('Pairwise interaction count')
dev.off()

