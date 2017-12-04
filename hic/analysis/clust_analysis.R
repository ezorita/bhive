suppressMessages(library(GenomicRanges))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

# Params
hic.binsize  = 5e3
expr.cutoff  = 0.4
enh.radius   = 2.5e3
prom.radius  = 2.5e3
hp.binsize   = 1e5
global.clust.names = list('0'='Act-','1'='Lamin','2'='Act+','3'='H3K9me','4'='H3K27me')

# Files
fig.dir      = 'figures/'
cmarks.files = 'chip_datasets/all/*.bed'
cluster.file = 'hi_c/3d_patterns.out'
hiv.datasets = 'hiv_datasets/*.integ'
hiv.expr     = 'expr_datasets/bhive_all_expr.txt'
gene.data    = 'annotations/Homo_Sapiens.GRCh37.75.gtf'
gene.expr    = 'expr_datasets/jurkat_rnaseq.tpm'
chrsz.data   = 'annotations/hg19_chromsizes.txt'
enh.data     = 'annotations/H3K27ac.01'
hotspot.data = 'hiv_datasets/jurkat_bhive_bushman_chronly_sorted.txt'
nup.damid    = 'chip_datasets/Nup153_IMR90.bed'

# Analysis setup
pattern_sizes_pie    = TRUE
pattern_category_pie = TRUE
intrachromosomal_chromatin_heatmaps = TRUE
interchromosomal_chromatin_heatmap  = TRUE
interchromosomal_hiv_heatmap        = TRUE


# Function handlers
setdiff_nw      <- function(...) { suppressWarnings(setdiff(...)) };
findOverlaps_nw <- function(...) { suppressWarnings(findOverlaps(...)) };


## File index
##
## 1. Read data.
## 1.1. Read chromatin marks.
## 1.2. Read 3D patterns.
## 1.3. Read HIV datasets.
## 1.4. Read annotations.
##
## 2. Data structures.
## 2.1. Compute genomic categories.
##
## 3. Signal distribution in 3D patterns.
## 3.1. Chromatin marks in 3D patterns.
## 3.2. HIV in 3D patterns.
##
## 


##
## 1. Read data.
##
cat('[1. Read data]\n')


##
## 1.1. Read chromatin marks.
##
cat('[1.1. Read chromatin marks]\n')
cm.gr    = GRanges()
cm.names = c()

dataset.id = 1
for (f in Sys.glob(cmarks.files)) {
    # Extract CM name from filename (text before '_')
    path     <- strsplit(f,'/')[[1]]
    fname    <- strsplit(path[length(path)],'.',fixed=T)[[1]][1]
    cm.name  <- strsplit(fname,'_',fixed=T)[[1]][1]
    cm.names <- c(cm.names,cm.name)
    # Read file (only chr[1-9XY][0-9]*)
    d       <- subset(read.table(f), grepl("^chr[0-9XY][0-9]*$",V1))
    d$V1    <- factor(d$V1, unique(d$V1))
    cm.gr   <- c(cm.gr, GRanges(seqnames=Rle(d$V1), ranges=IRanges(start=d$V2,end=d$V3),
                                dataset= rep(dataset.id,nrow(d))))
    # Verbose dataset info
    cat('[1.1] chromatin mark #',dataset.id,': ',cm.name,' (',nrow(d),' lines)\n',sep='')
    # Increase loop counter
    dataset.id <- dataset.id + 1
}


##
## 1.2. Read 3D patterns
##
cat('[1.2. Read 3D patterns]\n')
clust.gr = GRanges()

for (line in readLines(file(cluster.file))) {
    words      <- strsplit(line,split="\t")[[1]]
    # Parse columns.
    chr        <- words[1]
    local.clu  <- as.integer(words[2])
    global.clu <- as.integer(words[3])
    clu.idx    <- as.integer(unlist(strsplit(words[4], split=",")))
    # Assemble cluster GR object.
    clu.gr     <- GRanges(seqnames = chr, IRanges(start=clu.idx*hic.binsize, end=(clu.idx+1)*hic.binsize-1),
                          local.cluster  = rep(local.clu, length(clu.idx)),
                          global.cluster = rep(global.clust.names[[words[3]]], length(clu.idx)),
                          chrom          = rep(chr, length(clu.idx))
                          )
    clust.gr   <- c(clust.gr,clu.gr)
}
# Verbose clusters
cat('[1.2] ', length(unique(clust.gr$chrom)), ' chromosomes\n', sep='')
cat('[1.2] ', length(unique(clust.gr$local.cluster)), ' intrachromosomal 3D patterns\n', sep='')
cat('[1.2] ', length(unique(clust.gr$global.cluster)), ' global 3D patterns\n', sep='')


##
## 1.3. Read HIV datasets.
##
cat('[1.3. Read HIV datasets]\n')

## 1.3.1. Insertion sites.
cat('[1.3.1. HIV insertion sites]\n')
hiv.gr    = GRanges()
hiv.names = c()

dataset.id = 1
for (f in Sys.glob(hiv.datasets)) {
    # Extract dataset name from filename (text before '_')
    path      <- strsplit(f,'/')[[1]]
    fname     <- strsplit(path[length(path)],'.',fixed=T)[[1]][1]
    hiv.name  <- strsplit(fname,'_',fixed=T)[[1]][1]
    hiv.names <- c(hiv.names,hiv.name)
    # Read file (only chr[1-9XY][0-9]*)
    d         <- subset(read.table(f), grepl("^chr[0-9XY][0-9]*$",V1))
    d$V1      <- factor(d$V1, unique(d$V1))
    data.gr   <- GRanges(seqnames = Rle(d$V1), ranges = IRanges(start=d$V2, width=1), strand = Rle(d$V3),
                         dataset      = rep(dataset.id, nrow(d)),
                         dataset.name = rep(hiv.name, nrow(d))
                         )
    hiv.gr    <- c(hiv.gr,data.gr)
    # Verbose dataset info
    cat('[1.3.1] HIV dataset #',dataset.id,': ',hiv.name,' (',nrow(d),' insertion sites)\n',sep='')
    # Increase loop counter
    dataset.id <- dataset.id + 1
}

# 1.3.2. HIV expression (BHIVE)
cat('[1.3.2. HIV expression data (BHIVE)]\n')
bhive    <- read.table(hiv.expr, header=T)
bhive.gr <- GRanges(Rle(bhive$chr), IRanges(start=bhive$locus, width=1),
	strand=Rle(bhive$strand), brcd=bhive$brcd,
	nread=bhive$reads, mapq=bhive$mapq, rep=bhive$rep, expr=bhive$exprscore)
# Verbose
cat('[1.3.2] hiv insertions with expression info: ',sum(!is.na(bhive$exprscore)),'\n')

##
## 1.4. Read annotations.
##
cat('[1.4. Read annotations]\n')

# 1.4.1. Reference genome.
cat('[1.4.1. Reference genome]\n')
genome    <- read.table(chrsz.data, as.is=TRUE)
genome.gr <- GRanges(Rle(genome$V1), IRanges(start=1, width=genome$V2))
# Verbose genome
cat('[1.4.1] chromosomes: ',length(genome.gr),'\n',sep='')
cat('[1.4.1] coverage: ',sum(as.numeric(sum(coverage(genome.gr)))),' nt\n',sep='')

## 1.4.2. Read expression/gene data.
cat('[1.4.2. Gene expression info]\n')
genes     <- subset(read.table(gene.expr,header=T), type=='protein_coding')
genes$chr <- paste('chr',genes$chr,sep='')
cutoff    <- quantile(genes$tpm, expr.cutoff)
genact    <- subset(genes, tpm > cutoff)
gensil    <- subset(genes, tpm <= cutoff)
# Verbose RNA-seq info.
cat('[1.4.2] genes: ',nrow(genes),'\n',sep='')
cat('[1.4.2] expression cutoff: ',cutoff,' tpm\n',sep='')
cat('[1.4.2] active genes: ',nrow(genact),'\n',sep='')
cat('[1.4.2] silent genes: ',nrow(gensil),'\n',sep='')

# 1.4.3 Read enhancer data.
cat('[1.4.3. Enhancer data]\n')
enh           <- read.table(enh.data, comm="#")
enh.gr        <- reduce(GRanges(Rle(enh$V1), IRanges(start=enh$V2, enh$V3)))
start(enh.gr) <- start(enh.gr) - enh.radius
end(enh.gr)   <- end(enh.gr) + enh.radius
cat('[1.4.3] enhancer bed file: ',enh.data,'\n',sep='')
cat('[1.4.3] enriched intervals: ',length(enh.gr),'\n',sep='')
cat('[1.4.3] coverage: ',sum(as.numeric(sum(coverage(enh.gr)))),' nt\n',sep='')


##
## 2. Data structures.
##
cat('[2. Data structures]\n')


##
## 2.1. Compute genomic categories.
##
cat('[2.1. Genomic categories]\n')
genes.gr  <- GRanges(Rle(genes$chr), IRanges(start=genes$beg, end=genes$end),
                     strand=Rle(genes$strand), name=genes$gene_name, expr=genes$tpm)
actgen.gr <- GRanges(Rle(genact$chr), IRanges(start=genact$beg, end=genact$end),
                     strand=Rle(genact$strand), name=genact$gene_name, expr=genact$tpm)
silgen.gr <- GRanges(Rle(gensil$chr), IRanges(start=gensil$beg, end=gensil$end),
                     strand=Rle(gensil$strand), name=gensil$gene_name, expr=gensil$tpm)
actpro.gr <- promoters(actgen.gr, up=prom.radius, down=prom.radius)
acttss.gr <- resize(actgen.gr,1)
siltss.gr <- resize(silgen.gr,1)

## TODO: Why not new categories, like: gene+enhancer, promoter+enhacer, enhancer alone?
## TODO: Also for the promoter, would it make sense to split in internal, external? Or delete the category?

# Active genes prevail over silent ones.
silgen_na     <- setdiff_nw(silgen.gr, actgen.gr, ignore.strand = TRUE)
# Subtract promoter regions from genes.
actgen_np     <- setdiff_nw(actgen.gr, actpro.gr, ignore.strand = TRUE)
silgen_nanp   <- setdiff_nw(silgen_na, actpro.gr, ignore.strand = TRUE)
# Subtract enhancers from genes and promoters
actpro_ne     <- setdiff_nw(actpro.gr, enh.gr, ignore.strand = TRUE)
actgen_npne   <- setdiff_nw(actgen_np, enh.gr, ignore.strand = TRUE)
silgen_nanpne <- setdiff_nw(silgen_nanp, enh.gr, ignore.strand = TRUE)

# Create category ranges list
categ = list()
categ[['silgen']] = reduce(silgen_nanpne)
categ[['actgen']] = reduce(actgen_npne)
categ[['actpro']] = reduce(actpro_ne)
categ[['enh']]    = reduce(enh.gr)
categ[['int']]    = reduce(setdiff_nw(genome.gr,union(silgen_nanpne, union(actgen_npne, union(actpro_ne,enh.gr)) )))

# Verbose categories.
cat('[2.1] active genes: ',sum(as.numeric(sum(coverage(categ[['actgen']])))),' nt\n',sep='')
cat('[2.1] silent genes: ',sum(as.numeric(sum(coverage(categ[['silgen']])))),' nt\n',sep='')
cat('[2.1] active promoters:: ',sum(as.numeric(sum(coverage(categ[['actpro']])))),' nt\n',sep='')
cat('[2.1] enhancers: ',sum(as.numeric(sum(coverage(categ[['enh']])))),' nt\n',sep='')
cat('[2.1] intergenic: ',sum(as.numeric(sum(coverage(categ[['int']])))),' nt\n',sep='')

# Define category colors
categ.colors        <- c("#92C46DA0","#548B54", "#19485780", "#000000A0", "#C0955FE0")
names(categ.colors) <- c('actgen','actpro','enh','silgen','int')


##
## 2.2 3D patterns.
##
cat('[2.2. 3D patterns]\n')

# Cluster count and size.
clust.bins <- table(clust.gr$global.cluster)
clust.sz   <- clust.bins * hic.binsize
clust.cnt  <- length(clust.sz)

# Define 3D pattern colors
clust.colors               <- rev(brewer.pal(clust.cnt,"Blues"))
names(clust.colors)        <- c("Act+","Act-","H3K27me","H3K9me","Lamin")

# 3D pattern sizes pie chart.
if (pattern_sizes_pie) {
    figure.fn <- paste(fig.dir,'3dpattern_sizes.pdf',sep='')
    cat('[2.2] cluster sizes -> ',figure.fn,'\n',sep='')

    # Pie chart, cluster sizes.
    pdf(figure.fn,useDingbats=F)
    labels = paste('Cluster ',names(clust.sz),'\n',clust.sz/1e6,' Mbp')
    pie(clust.bins,labels=labels,col=clust.colors[names(clust.bins)],main='3D pattern sizes')
    dev.off()
}


##
## 2.3 Category background of 3D patterns
##
cat('[2.3 Category background of 3D patterns]\n')

# Compute 3D pattern / category intersect
clust.categ    <- data.frame(clust=c(),cat=c(),span=c())
clust.categ.gr <- GRanges()
for (i in unique(clust.gr$global.cluster)) {
    clust.df <- data.frame(clust=c(),cat=c(),span=c())
    # iterate over pattern categories
    for (idx in names(categ)) {
        gc             <- clust.gr[clust.gr$global.cluster == i]
        # Intersect and annotate cluster and background category.
        clu.cat        <- intersect(gc,categ[[idx]],ignore.strand=T)
        clu.cat$cat    <- idx
        clu.cat$clust  <- i
        clust.categ.gr <- c(clust.categ.gr, clu.cat)
        # Store intersect span.
        categ.df       <- data.frame(clust=toString(i),cat=idx,span=sum(as.numeric(sum(coverage(clu.cat)))))
        categ.df$cat   <- levels(droplevels(categ.df$cat))
        clust.df       <- rbind(clust.df,categ.df)
    }
    clust.categ <- rbind(clust.categ,clust.df)
}

# Plot background categories in pie charts.
if (pattern_category_pie) {
    figure.fn = paste(fig.dir,'3dpattern_background_category.pdf',sep='')
    cat('[2.3] 3D pattern background category pie charts -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=20,height=15)
    par(mfrow=c(2,ceiling(clust.cnt/2)))
    for (i in unique(clust.gr$global.cluster)) {
        clust.df <- clust.categ[clust.categ$clust == i,]
        pie(clust.df$span, labels=paste(clust.df$cat,'\n',round(clust.df$span/1e4)/1e2,' Mbp',sep=''),main=paste('cluster',i),col=categ.colors[clust.df$cat])
    }
    dev.off()
}


##
## 3. Signal distribution in 3D patterns.
##
cat('[3. Signal in 3D patterns]\n')


##
## 3.1. Chromatin marks in 3D patterns.
##
cat('[3.1. Chromatin marks in 3D patterns]\n')

# Preallocate matrix
mat = list();
for (chr in levels(seqnames(clust.gr)))
    mat[[chr]] <- matrix(0,ncol=max(clust.gr[seqnames(clust.gr) == chr]$local.cluster),nrow=length(names))
gwmat <- matrix(0,ncol=length(unique(clust.gr$global.cluster)),nrow=length(names))

# Heatmap colors
hmcols <- colorRampPalette(c("blue","white","red"))(32)

# Local intrachromosomal clusters
normat <- list();
for (chr in levels(seqnames(clust.gr))) {
    # Select chromosome
    gc = clust.gr[seqnames(clust.gr) == chr]
    gd = cm.gr[seqnames(cm.gr) == chr]

    # Overlap Genomic Ranges
    ov = findOverlaps_nw(gc,gd,type='any')
    
    # Fill matrices
    mat[[chr]] = as.matrix(table(data.frame(dataset=gd[ov@to]$dataset,cluster=gc[ov@from]$local.cluster)))
    clsz = table(gc$local.cluster)
    normat[[chr]] = mat[[chr]]/matrix(clsz,nrow=dim(mat[[chr]])[1],ncol=length(clsz),byrow=TRUE)/rowSums(mat[[chr]])

    # Intrachromosomal cluster figure
    if (intrachromosomal_chromatin_heatmaps) {
        figure.fn = paste(fig.dir,'heatmap_3dpattern_chrommarks_',chr,'.pdf',sep='')
        cat('[3.1] intrachromosomal ',chr,' -> ',figure.fn,'\n',sep='')
        pdf(figure.fn,useDingbats=F,width=12,height=8)
        # Add 'trace = "none"' for presentation figure.
        heatmap.2(normat[[chr]],dendrogram='both',symbreaks=F,labRow=cm.names,scale='row',col=hmcols,margin=c(4,8),
                  symkey=F,key.title="Color key",xlab='3D patterns',
                  main=paste('Jurkat',chr,'3D patterns vs Chromatin marks'))
        dev.off()
    }
}

# Global cluster vs chromatin marks matrix
clu.cm.ov <- findOverlaps_nw(clust.gr,cm.gr,type='any')
clu.cm.df <- data.frame(dataset=cm.gr[clu.cm.ov@to]$dataset,cluster=clust.gr[clu.cm.ov@from]$global.cluster)
clu.cm.gw <- as.matrix(table(clu.cm.df))

# Global cluster vs chromatin marks heatmap
if (interchromosomal_chromatin_heatmap) {
    # Normalize matrix by rowSums
    clu.cm.mat  <- clu.cm.gw/matrix(clust.bins,nrow=dim(clu.cm.gw)[1],ncol=length(clust.bins),byrow=TRUE)/rowSums(clu.cm.gw)
    # Create heatmap
    clu.cm.hmap <- heatmap.2(clu.cm.mat,dendrogram='both',symbreaks=F,labRow=cm.names,scale='row',col=hmcols,
                                margin=c(4,8),symkey=F,key.title="Color key",xlab='3D pattern',
                                main='Jurkat GW 3D patterns vs Chromatin marks')
    clu.cm.lab  <- colnames(clu.cm.mat)[clu.cm.hmap$colInd]
    # Save figure
    figure.fn <- paste(fig.dir,'heatmap_3dpattern_chrommarks.pdf',sep='')
    cat('[3.1] heatmap: chromatin marks vs 3D patterns -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=12,height=15)
    # Add 'trace = "none"' for presentation figure
    heatmap.2(clu.cm.mat, dendrogram='both', symbreaks=FALSE, labRow=cm.names,scale='row', col=hmcols,margin=c(4,8),
              symkey=FALSE, key.title="Color key", xlab='3D pattern', main='Jurkat GW 3D patterns vs Chromatin marks',
              add.expr = text(x = seq_along(clu.cm.lab), y = 0.1, cex=1.5, srt =0, labels=clu.cm.lab, xpd=TRUE),
              labCol="")
    dev.off()
}


##
## 3.2. HIV in 3D patterns.
##
cat('[3.2. HIV insertion sites in 3D patterns]\n')

# Heatmap of HIV insertion sites in 3D patterns
clu.hiv.ov <- findOverlaps_nw(clust.gr,hiv.gr,type='any',ignore.strand=T)
clu.hiv.df <- data.frame(dataset=hiv.gr[clu.hiv.ov@to]$dataset.name, cluster=clust.gr[clu.hiv.ov@from]$global.cluster)
clu.hiv.gw <- as.matrix(table(clu.hiv.df))

if (interchromosomal_hiv_heatmap) {
    # Normalize matrix by rowSums
    clu.hiv.mat  <- clu.hiv.gw/matrix(clust.bins,nrow=dim(clu.hiv.gw)[1],ncol=length(clust.bins),byrow=TRUE)/rowSums(clu.hiv.gw)
    # Create heatmap
    clu.hiv.hmap <- heatmap.2(clu.hiv.mat,dendrogram='column',symbreaks=F,labRow=hiv.names,scale='row',
                                 col=hmcols,margin=c(10,20),symkey=F,key.title="Color key",xlab='3D pattern',
                                 main='Jurkat GW 3D patterns vs HIV insertion sites')
    clu.hiv.lab  <- colnames(clu.hiv.mat)[clu.hiv.hmap$colInd]
    # Save figure
    figure.fn <- paste(fig.dir,'heatmap_3dpattern_hiv.pdf',sep='')
    cat('[3.2] heatmap: HIV insertion sites vs 3D patterns -> ',figure.fn,'\n',sep='')
    pdf(figure.fn,useDingbats=F,width=16,height=5)
    # Add 'trace = "none"' for presentation figure
    heatmap.2(clu.hiv.mat, dendrogram='column', symbreaks=FALSE, labRow=hiv.names,scale='row', col=hmcols,
              margin=c(10,20), key=FALSE, xlab='3D pattern', main='Jurkat GW 3D patterns vs HIV insertion sites',
              add.expr = text(x = seq_along(clu.hiv.lab), y = 0.1, cex=2, srt =0, labels=clu.hiv.lab,
                              xpd=TRUE),
              labCol="")
    dev.off()
}


##
## 4. HIV integration pattern
##
cat('[4. HIV integration pattern]\n');

# Compute HIV frequencies in cluster/categories.
clu.cat.hiv.ov    <- findOverlaps_nw(hiv.gr, clust.categ.gr, ignore.strand=TRUE)
clu.cat.hiv       <- as.data.frame(table(clust.categ.gr[clu.cat.hiv.ov@to]@elementMetadata))
clu.cat.hiv       <- merge(clust.categ,clu.cat.hiv)
# Sort categories and clusters
clu.cat.hiv$cat   <- factor(clu.cat.hiv$cat, c("actgen","actpro","enh","silgen","int"))
clu.cat.hiv$clust <- factor(clu.cat.hiv$clust, c("Act+","Act-","H3K27me","H3K9me","Lamin"))
clu.cat.hiv       <- clu.cat.hiv[order(clu.cat.hiv$clu, clu.cat.hiv$cat),]


#
# 4.1. Integration density in Categories.
#
cat('[4.1. Integration density bias in Categories]\n')

# Plot spie chart of categories.
cat.hiv = aggregate(cbind(span, Freq) ~ cat, data=clu.cat.hiv, FUN=sum)

# GSpie shapes
cat.hiv$xmax <- cumsum(cat.hiv$span)
cat.hiv$xmin <- c(0, cat.hiv[1:(nrow(cat.hiv)-1),]$xmax)
cat.hiv$ymin <- 0
cat.hiv$ymax <- sqrt( cat.hiv$Freq/sum(cat.hiv$Freq) / (cat.hiv$span / sum(cat.hiv$span)) )
# Plot spie-chart using ggplot2
figure.fn <- paste(fig.dir,'hiv_bias_category.pdf',sep='')
cat('[4.1] HIV integration bias in categories -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=7,height=7)
ggplot(cat.hiv) +
    geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
    geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
    geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color=categ.colors[cat.hiv$cat]) +
    scale_fill_manual(name='Category',values=categ.colors) +
    coord_polar('x') +
    theme_minimal() +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    ggtitle('HIV integration bias in Categories')
dev.off()


#
# 4.2. Integration density in 3D patterns.
#
cat('[4.2. Integration density bias in 3D Patterns]\n')

# Plot spie chart of categories.
clu.hiv = aggregate(cbind(span, Freq) ~ clust, data=clu.cat.hiv, FUN=sum)

# GSpie shapes
clu.hiv$xmax <- cumsum(clu.hiv$span)
clu.hiv$xmin <- c(0, clu.hiv[1:(nrow(clu.hiv)-1),]$xmax)
clu.hiv$ymin <- 0
clu.hiv$ymax <- sqrt( clu.hiv$Freq/sum(clu.hiv$Freq) / (clu.hiv$span / sum(clu.hiv$span)) )
# Plot spie-chart using ggplot2
figure.fn <- paste(fig.dir,'hiv_bias_3dpattern.pdf',sep='')
cat('[4.2] HIV integration bias in 3D Patterns -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=7,height=7)
ggplot(clu.hiv) +
    geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
    geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
    geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color='black') +
    scale_fill_manual(name='3D Pattern',values=clust.colors) +
    coord_polar('x') +
    theme_minimal() +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    ggtitle('HIV integration bias in 3D Patterns')
dev.off()


#
# 4.3. Local integration density in categories (grouped and normalized by 3D pattern)
#
cat('[4.3. Integration density bias in categories (grouped and normalized by 3D pattern)]\n')

# Compute signal/background distributions and spie shapes.
spie.list = list()
for (i in unique(clust.gr$global.cluster)) {
    cat.hiv <- clu.cat.hiv[clu.cat.hiv$clust == i,]
    # GSpie shapes
    cat.hiv$xmax <- cumsum(cat.hiv$span)
    cat.hiv$xmin <- c(0, cat.hiv[1:(nrow(cat.hiv)-1),]$xmax)
    cat.hiv$ymin <- 0
    cat.hiv$ymax <- sqrt( cat.hiv$Freq/sum(cat.hiv$Freq) / (cat.hiv$span / sum(cat.hiv$span)) )
    # Plot spie-chart using ggplot2
    spie.plt <- ggplot(cat.hiv) +
        geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
        geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
        geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color=categ.colors[cat.hiv$cat]) +
        scale_fill_manual(name='Category',values=categ.colors) +
        coord_polar('x') +
        theme_minimal() +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
              axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        ggtitle(i)
    # Append to plot list
    spie.list[[i]] = spie.plt
}

# Plot HIV category biases in spie-charts.
figure.fn <- paste(fig.dir,'hiv_bias_categories_in_3dpatterns.pdf',sep='')
cat('[4.3] HIV local integration bias (category) in 3D patterns -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=20,height=15)
do.call("grid.arrange", c(spie.list, ncol=floor(sqrt(length(spie.list)))))
dev.off()


#
# 4.4. Local integration density in 3D patterns (grouped and normalized by Category)
#
cat('[4.4. Integration density bias in 3D pattern (grouped and normalized by Category)]\n')

# Compute signal/background distributions and spie shapes.
spie.list = list()
for (i in names(categ)) {
    clu.hiv <- clu.cat.hiv[clu.cat.hiv$cat == i,]
    # GSpie shapes
    clu.hiv$xmax <- cumsum(clu.hiv$span)
    clu.hiv$xmin <- c(0, clu.hiv[1:(nrow(clu.hiv)-1),]$xmax)
    clu.hiv$ymin <- 0
    clu.hiv$ymax <- sqrt( clu.hiv$Freq/sum(clu.hiv$Freq) / (clu.hiv$span / sum(clu.hiv$span)) )
    # Plot spie-chart using ggplot2
    spie.plt <- ggplot(clu.hiv) +
        geom_hline(yintercept=sqrt(c(1)),alpha=0.5) +
        geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") +
        geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color='black') +
        scale_fill_manual(name='3D pattern',values=clust.colors) +
        coord_polar('x') +
        theme_minimal() +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
              axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        ggtitle(i)
    # Append to plot list
    spie.list[[i]] = spie.plt
}

# Plot HIV category biases in spie-charts.
figure.fn <- paste(fig.dir,'hiv_bias_3dpatterns_in_categories.pdf',sep='')
cat('[4.4] HIV local integration bias (3D pattern) in categories -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=20,height=15)
do.call("grid.arrange", c(spie.list, ncol=floor(sqrt(length(spie.list)))))
dev.off()


#
# 4.5. HIV integration density (grouped by 3D pattern or category)
#
cat('[4.5. HIV integration density in 3D patterns and categories]\n')

# Compute integration density (per Mb) in cluster+category.
clu.cat.hiv$integ.mb <- clu.cat.hiv$Freq/clu.cat.hiv$span*1e6

# Barplot of 3d patterns grouped by category.
figure.fn <- paste(fig.dir,'hiv_density_category_3dpattern.pdf',sep='')
cat('[4.5] HIV integration density (groups: category) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=8,height=5)
ggplot(data=clu.cat.hiv,aes(fill=clust,x=cat,y=integ.mb)) +
    geom_bar(stat='identity',position='dodge',colour='black') +
    scale_fill_manual(name='3D pattern',values=c(clust.colors)) +
    xlab('Category') +
    ylab('HIV integrations per Mb')
dev.off()

# Barplot of categories grouped by 3d pattern.
figure.fn <- paste(fig.dir,'hiv_density_3dpattern_category.pdf',sep='')
cat('[4.5] HIV integration density (groups: 3D pattern) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=8,height=5)
ggplot(data=clu.cat.hiv,aes(fill=cat,x=clust,y=integ.mb)) +
    geom_bar(stat='identity',position='dodge',colour='black') +
    scale_fill_manual(name='Category',values=c(categ.colors)) +
    xlab('3D patterns') +
    ylab('HIV integrations per Mb')
dev.off()



##
## 5. Gene expression in 3D patterns
##
cat('\n[5. Gene expression in 3D patterns]\n')


#
# 5.1 Endogenous expression in 3D patterns
#
cat('[5.1. Endogenous expression in 3D patterns]\n')
genes.expr.ov  <- findOverlaps_nw(clust.gr, acttss.gr,ignore.strand=T)
genes.expr.clu <- data.frame(clust = clust.gr[genes.expr.ov@from]$global.cluster,
                             expr  = acttss.gr[genes.expr.ov@to]$expr)
# Sort 3D Pattern factors
genes.expr.clu$clust <- factor(genes.expr.clu$clust,c("Act+","Act-","H3K9me","H3K27me","Lamin"))

# Barplot of endogenous expression in 3D patterns.
figure.fn <- paste(fig.dir,'expr_endog_3dpattern.pdf',sep='')
cat('[5.1] Endogenous expression distribution in 3D patterns (boxplot) -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=7, height=7)
ggplot(genes.expr.clu, aes(x = clust, y = expr, fill = clust)) +
    geom_boxplot(outlier.size=0.3) +
    scale_y_continuous(trans='log10') +
    scale_fill_manual(name='',values = clust.colors) +
    theme_minimal() +
    xlab('') +
    ylab('Expression of Jurkat active genes [tpm]')
dev.off()


#
# 5.2 HIV expression in 3D patterns
#
cat('[5.2. HIV expression in 3D patterns]\n')

# Compute clusters/categories of BHIVE integrations.
clu.cat.bhive.ov <- findOverlaps_nw(bhive.gr, clust.categ.gr, ignore.strand=TRUE)
clu.cat.bhive    <- data.frame(clust = clust.categ.gr[clu.cat.bhive.ov@to]$clust,
                               cat   = clust.categ.gr[clu.cat.bhive.ov@to]$cat,
                               expr   = bhive.gr[clu.cat.bhive.ov@from]$expr)
# Sort categories and clusters
clu.cat.bhive$cat   <- factor(clu.cat.bhive$cat, c("actgen","actpro","enh","silgen","int"))
clu.cat.bhive$clust <- factor(clu.cat.bhive$clust, c("Act+","Act-","H3K27me","H3K9me","Lamin"))
clu.cat.bhive       <- clu.cat.bhive[order(clu.cat.bhive$clu, clu.cat.bhive$cat),]

# Expression distribution by clusters.
p1 <- ggplot(clu.cat.bhive, aes(x = cat, y = expr, fill = cat)) +
    geom_boxplot(outlier.size=0.3) +
    scale_fill_manual(name='',values = categ.colors) +
    coord_cartesian(ylim=c(0,3)) +
    theme_minimal() +
    xlab('') +
    ylab('Normalized HIV expression (log10)')

# Expression distribution by categories.
p2 <- ggplot(clu.cat.bhive, aes(x = clust, y = expr, fill = clust)) +
    geom_boxplot(outlier.size=0.3) +
    scale_fill_manual(name='',values = clust.colors) +
    coord_cartesian(ylim=c(0,3)) +
    theme_minimal() +
    xlab('') +
    ylab('')

# Expression by cluster and category (Grouped by cluster)
p3 <- ggplot(clu.cat.bhive, aes(x = clust, y = expr, fill = cat)) +
    geom_boxplot(outlier.size=0.3) +
    scale_fill_manual(name='Category',values = categ.colors) +
    coord_cartesian(ylim=c(0,3)) +
    theme_minimal() +
    xlab('3D pattern') +
    ylab('Normalized HIV expression (log10)')

# Expression by cluster and category (Grouped by category)
p4 <- ggplot(clu.cat.bhive, aes(x = cat, y = expr, fill = clust)) +
    geom_boxplot(outlier.size=0.3) +
    scale_fill_manual(name='3D pattern',values = clust.colors) +
    coord_cartesian(ylim=c(0,3)) +
    theme_minimal() +
    xlab('Category') +
    ylab('')

figure.fn <- paste(fig.dir,'expr_bhive.pdf',sep='')
cat('[5.2] HIV expression in 3D patterns (full figure) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=20,height=15)
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()

##
## END OF REVISED FILE ##
##


categ.expr = data.frame(clust=c(),cat=c(),expr=c())
clust.full = data.frame(clust=c(),cat=c(),span=c(),cat.ratio=c(),integ=c(),local.expected=c(),global.expected=c(),integ.mb=c(),local.enrich=c(),global.enrich=c())
# Endognous expression.
endog.expr = data.frame(clust=c(),expr=c())
gactgen.tss = promoters(actgen,up=1,down=1)
global.ratio = nrow(hiv)/sum(clust.categ$span)
# Integration bias spie-chart (2/4).


    clust.df = clust.categ[clust.categ$clust == toString(i),]
    clust.df$cat.ratio = clust.df$span/sum(clust.df$span)
    lat.df = clust.df
    clust.df$integ = 0
    gc = clust.gr[clust.gr$global.cluster == toString(i)]
    for (idx in names(categ)) {
        gcat = intersect(gc,categ[[idx]],ignore.strand=T)
        # Compute HIV expression in Cluster categories
        ov = suppressWarnings(findOverlaps(gcat,ghiv,ignore.strand=T))
        expr.df = data.frame(clust=toString(i),cat=idx,expr=ghiv[ov@to]$expr)
        expr.df$cat = levels(droplevels(expr.df$cat))
        categ.expr = rbind(categ.expr, expr.df)
        clust.df[clust.df$cat == idx,]$integ = length(unique(ov@to))
    }
    # Global clusters all integ datasets
    clust.df$local.expected = clust.df$cat.ratio*sum(clust.df$integ)
    clust.df$global.expected = global.ratio*clust.df$span
    clust.df$integ.mb = clust.df$integ / clust.df$span * 1e6
    clust.df$local.enrich = clust.df$integ / clust.df$local.expected
    clust.df$global.enrich = clust.df$integ / clust.df$global.expected
    clust.full = rbind(clust.full,clust.df)
    # Integration bias spie-chart(3/4).
    spiechart(expected=clust.df$local.expected,observed=clust.df$integ, col=categ.colors[clust.df$cat], main=paste('cluster',i),labels=paste(clust.df$cat,'(',clust.df$integ,'/',clust.df$expected,')',sep=''))
    # Compute Endogenous expression in Cluster categories
    ov.end = suppressWarnings(findOverlaps(gc,gactgen.tss,ignore.strand=T))
    if (length(ov.end) > 0) {
        endog.expr.df = data.frame(clust=toString(i),expr=gactgen.tss[ov.end@to]$expr)
        endog.expr = rbind(endog.expr,endog.expr.df)
    }
}
# Integration bias spie-chart(4/4).
dev.off()





#   - Load Nup data.
nup.data = read.table(nup.damid)
gnup = GRanges(seqnames=Rle(nup.data$V1),ranges=IRanges(start=nup.data$V2,end=nup.data$V3))

# NAD categories
nad.categ = list()
nad.categ[['silgen-nad']] = suppressWarnings(intersect(categ[['silgen']],gnup))
nad.categ[['actgen-nad']] = suppressWarnings(intersect(categ[['actgen']],gnup))
nad.categ[['actpro-nad']] = suppressWarnings(intersect(categ[['actpro']],gnup))
nad.categ[['enh-nad']] = suppressWarnings(intersect(categ[['enh']],gnup))
nad.categ[['int-nad']] = suppressWarnings(intersect(categ[['int']],gnup))
# Non-NAD categories
nad.categ[['silgen']] = setdiff_nw(categ[['silgen']],gnup)
nad.categ[['actgen']] = setdiff_nw(categ[['actgen']],gnup)
nad.categ[['actpro']] = setdiff_nw(categ[['actpro']],gnup)
nad.categ[['enh']] = setdiff_nw(categ[['enh']],gnup)
nad.categ[['int']] = setdiff_nw(categ[['int']],gnup)

# SuperSpyChart Nup153
#   - Add category info to clusters.
clust.gr$cat = ''
#   * The order of the loop is important to correctly resolve bins falling in category edges.
for (idx in c('int','silgen','actgen','actpro','enh')) {
    ov = suppressWarnings(findOverlaps(clust.gr,categ[[idx]]))
    clust.gr[ov@from]$cat = idx
}

nup.all = data.frame(clust=c(),cat=c(),bg.nt=c(),nup.nt=c())
for (i in unique(clust.gr$global.cluster)) {
    for (idx in names(categ)) {
        gc = clust.gr[clust.gr$global.cluster == i & clust.gr$cat == idx]
        nup.df = data.frame(clust=toString(i),cat=idx,bg.nt=sum(as.numeric(sum(coverage(gc)))),nup.nt=sum(as.numeric(sum(coverage(intersect(gc,gnup,ignore.strand=T))))))
        nup.df$cat = levels(droplevels(nup.df$cat))
        nup.all = rbind(nup.all,nup.df)
    }
}

#   Spie on category.
nup.cat = merge(aggregate(bg.nt ~ cat,data=nup.all,FUN=sum),aggregate(nup.nt ~ cat,data=nup.all,FUN=sum),by='cat')
nup.cat$cat = factor(nup.cat$cat,c("actgen","actpro","enh","silgen","int"))
nup.cat = nup.cat[order(nup.cat$cat),]
nup.cat$xmax = cumsum(nup.cat$bg.nt)
nup.cat$xmin = c(0,nup.cat[1:(nrow(nup.cat)-1),]$xmax)
nup.cat$ymin = 0
nup.cat$ymax = sqrt(nup.cat$nup.nt/sum(nup.cat$nup.nt)/(nup.cat$bg.nt/sum(nup.cat$bg.nt)))
pdf(paste(fig.dir,'nup_category_signal_spie.pdf',sep=''),useDingbats=F,width=8,height=8)
ggplot(nup.cat) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color=categ.colors[nup.cat$cat]) + scale_fill_manual(name='Category',values=categ.colors)+ coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Nup153 enrichment')
dev.off()
#   Spie on 3D pattern.
nup.clu = merge(aggregate(bg.nt ~ clust,data=nup.all,FUN=sum),aggregate(nup.nt ~ clust,data=nup.all,FUN=sum),by='clust')
nup.clu$clust = factor(nup.clu$clust,c("Act+","Act-","H3K27me","H3K9me","Lamin"))
nup.clu = nup.clu[order(nup.clu$clust),]
nup.clu$xmax = cumsum(nup.clu$bg.nt)
nup.clu$xmin = 0
nup.clu[2:nrow(nup.clu),]$xmin = nup.clu[1:(nrow(nup.clu)-1),]$xmax
nup.clu$ymin = 0
nup.clu$ymax = sqrt(nup.clu$nup.nt/sum(nup.clu$nup.nt)/(nup.clu$bg.nt/sum(nup.clu$bg.nt)))
pdf(paste(fig.dir,'nup_clust_signal_spie.pdf',sep=''),useDingbats=F,width=8,height=8)
ggplot(nup.clu) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color="#555555") + scale_fill_manual(name='3D pattern',values=clust.colors)+ coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Nup153 enrichment')
dev.off()

# Crazy barplot
nup.all$clust = factor(nup.all$clust,c("Act+","Act-","H3K27me","H3K9me","Lamin"))
nup.all$cat = factor(nup.all$cat,c("actgen","actpro","enh","silgen","int"))
nup.all = nup.all[order(nup.all$cat,nup.all$clust),]
nup.all$ratio.bg = nup.all$bg.nt / sum(nup.all$bg.nt)
nup.all$ratio.nup = nup.all$nup.nt / sum(nup.all$nup.nt)
nup.all$xmax = cumsum(nup.all$bg.nt)
nup.all$xmin = 0
nup.all[2:nrow(nup.all),]$xmin = nup.all[1:(nrow(nup.all)-1),]$xmax
nup.all$ymax = sqrt(nup.all$ratio.nup/nup.all$ratio.bg)
nup.all$ymin = 0

nup.sq = nup.cat[,c('cat','xmin','xmax','ymin','ymax')]
nup.sq$clust = nup.sq$cat
for (c in unique(nup.sq$cat)) {
    nup.all[nup.all$cat == c,'ymin'] = nup.sq[nup.sq$cat == c,'ymax']
}
nup.sq = rbind(nup.sq,nup.all[,c('cat','xmin','xmax','ymin','ymax','clust')])
clust.colors.gray = rep("#555555",5)
names(clust.colors.gray) = names(clust.colors)
pdf(paste(fig.dir,'nup_clust_cat_signal_superspie.pdf',sep=''),useDingbats=F,width=12,height=12)
ggplot(nup.sq) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax,color=clust),alpha=0.5) + scale_fill_manual(name='',values=c(clust.colors,categ.colors)) + scale_color_manual(name='',values=c(clust.colors.gray,categ.colors))+ coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Nup153 enrichment by Category and 3D pattern')
dev.off()

# Nup Barplot all. (it does not reflect well the reality beacuse the bars may represent tiny or huge regions and are all the same width).
nup.all$cov = nup.all$nup.nt/nup.all$bg.nt
ggplot(nup.all,aes(fill=clust,x=cat,y=cov)) + geom_bar(stat='identity',position="dodge",color="#555555") + scale_fill_manual(values=clust.colors)

# Integration bias spie-chart (1/4).
source('spie.R')
pdf(paste(fig.dir,'global_cluster_hiv_bias.pdf',sep=''),useDingbats=F,width=20,height=15)
categ.expr = data.frame(clust=c(),cat=c(),expr=c())
clust.full = data.frame(clust=c(),cat=c(),span=c(),cat.ratio=c(),integ=c(),local.expected=c(),global.expected=c(),integ.mb=c(),local.enrich=c(),global.enrich=c())
# Endognous expression.
endog.expr = data.frame(clust=c(),expr=c())
gactgen.tss = promoters(actgen,up=1,down=1)
global.ratio = nrow(hiv)/sum(clust.categ$span)
# Integration bias spie-chart (2/4).
par(mfrow=c(2,ceiling(clust.cnt/2)))
for (i in unique(clust.gr$global.cluster)) {
    clust.df = clust.categ[clust.categ$clust == toString(i),]
    clust.df$cat.ratio = clust.df$span/sum(clust.df$span)
    clust.df$integ = 0
    gc = clust.gr[clust.gr$global.cluster == toString(i)]
    for (idx in names(categ)) {
        gcat = intersect(gc,categ[[idx]],ignore.strand=T)
        # Compute HIV expression in Cluster categories
        ov = suppressWarnings(findOverlaps(gcat,ghiv,ignore.strand=T))
        expr.df = data.frame(clust=toString(i),cat=idx,expr=ghiv[ov@to]$expr)
        expr.df$cat = levels(droplevels(expr.df$cat))
        categ.expr = rbind(categ.expr, expr.df)
        clust.df[clust.df$cat == idx,]$integ = length(unique(ov@to))
    }
    # Global clusters all integ datasets
    clust.df$local.expected = clust.df$cat.ratio*sum(clust.df$integ)
    clust.df$global.expected = global.ratio*clust.df$span
    clust.df$integ.mb = clust.df$integ / clust.df$span * 1e6
    clust.df$local.enrich = clust.df$integ / clust.df$local.expected
    clust.df$global.enrich = clust.df$integ / clust.df$global.expected
    clust.full = rbind(clust.full,clust.df)
    # Integration bias spie-chart(3/4).
    spiechart(expected=clust.df$local.expected,observed=clust.df$integ, col=categ.colors[clust.df$cat], main=paste('cluster',i),labels=paste(clust.df$cat,'(',clust.df$integ,'/',clust.df$expected,')',sep=''))
    # Compute Endogenous expression in Cluster categories
    ov.end = suppressWarnings(findOverlaps(gc,gactgen.tss,ignore.strand=T))
    if (length(ov.end) > 0) {
        endog.expr.df = data.frame(clust=toString(i),expr=gactgen.tss[ov.end@to]$expr)
        endog.expr = rbind(endog.expr,endog.expr.df)
    }
}
# Integration bias spie-chart(4/4).
dev.off()

# NAD integration bias vs category.
nad.bias.cat = data.frame(cat=c(),span=c(),cat.ratio=c(),integ=c(),integ.expected=c(),integ.mb=c())
total.span = sum(clust.categ$span)
for (i in names(nad.categ)) {
    gcat = nad.categ[[i]]
    ov.cnt = suppressWarnings(countOverlaps(gcat,ghiv,ignore.strand=T))
    span = sum(as.numeric(sum(coverage(gcat))))
    intexpect = nrow(hiv) * span/total.span
    integ = sum(ov.cnt)
    nad.bias.cat = rbind(nad.bias.cat, data.frame(cat=i,span=span,cat.ratio=span/total.span,integ=integ,integ.expected=intexpect, integ.mb=integ*1e6/span))
}
# Prepare data (xmin,xmax,ymin,ymax) for spie chart.
nad.bias.cat$cat = factor(nad.bias.cat$cat,c("actgen/NAD","actgen","actpro/NAD","actpro","enh/NAD","enh","silgen/NAD","silgen","int/NAD","int"))
nad.bias.cat = nad.bias.cat[order(nad.bias.cat$cat),]
nad.bias.cat$xmax = cumsum(nad.bias.cat$span)
nad.bias.cat$xmin = 0
nad.bias.cat[2:nrow(nad.bias.cat),]$xmin = nad.bias.cat[1:(nrow(nad.bias.cat)-1),]$xmax
nad.bias.cat$ymin = 0
nad.bias.cat$ymax = sqrt(nad.bias.cat$integ/nad.bias.cat$integ.expected)

# NAD spie colors (category)
nad.cat.colors = rep(c("#92C46DA0","#548B54","#000000A0", "#19485780", "#C0955FE0"),2)
names(nad.cat.colors) = c('actgen','actpro','silgen','enh','int','actgen/NAD','actpro/NAD','silgen/NAD','enh/NAD','int/NAD')
nad.cat.lines = c(rep('blank',5),rep('solid',5))
names(nad.cat.lines) = c('actgen','actpro','silgen','enh','int','actgen/NAD','actpro/NAD','silgen/NAD','enh/NAD','int/NAD')

# Generate spie with polar coordinates.
pdf(paste(fig.dir,'nup_category_integ_spie.pdf',sep=''),useDingbats=F,width=7,height=7)
ggplot(nad.bias.cat) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5) + coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Integration enrichment in NAD (Nup153 Associated Domain) / Categories') + scale_linetype_manual(name='Category',values=nad.cat.lines) + scale_fill_manual(name='Category',values=nad.cat.colors)
dev.off()

# NAD integration bias vs 3d pattern.
nad.bias.clust = data.frame(clust=c(),span=c(),clust.ratio=c(),integ=c(),integ.expected=c(),integ.mb=c())
for (i in unique(clust.gr$global.cluster)) {
    gclu = clust.gr[clust.gr$global.cluster == toString(i),]
    gclu.nad = suppressWarnings(intersect(gclu,gnup,ignore.strand=T))
    gclu.non = setdiff_nw(gclu,gnup,ignore.strand=T)
    nad.cnt = sum(suppressWarnings(countOverlaps(gclu.nad,ghiv,ignore.strand=T)))
    non.cnt = sum(suppressWarnings(countOverlaps(gclu.non,ghiv,ignore.strand=T)))
    nad.span = sum(as.numeric(sum(coverage(gclu.nad))))
    non.span = sum(as.numeric(sum(coverage(gclu.non))))
    nad.intexpect = nrow(hiv) * nad.span/total.span
    non.intexpect = nrow(hiv) * non.span/total.span
    nad.bias.clust = rbind(nad.bias.clust, data.frame(clust=paste(i,"/NAD",sep=""),span=nad.span,clust.ratio=nad.span/total.span,integ=nad.cnt,integ.expected=nad.intexpect, integ.mb=nad.cnt*1e6/span))
    nad.bias.clust = rbind(nad.bias.clust, data.frame(clust=i,span=non.span,clust.ratio=non.span/total.span,integ=non.cnt,integ.expected=non.intexpect, integ.mb=non.cnt*1e6/span))
}
# Prepare data (xmin,xmax,ymin,ymax) for spie chart.
nad.bias.clust$clust = factor(nad.bias.clust$clust,c("Act+/NAD","Act+","Act-/NAD","Act-","H3K27me/NAD","H3K27me","H3K9me/NAD","H3K9me","Lamin/NAD","Lamin"))
nad.bias.clust = nad.bias.clust[order(nad.bias.clust$clust),]
nad.bias.clust$xmax = cumsum(nad.bias.clust$span)
nad.bias.clust$xmin = 0
nad.bias.clust[2:nrow(nad.bias.clust),]$xmin = nad.bias.clust[1:(nrow(nad.bias.clust)-1),]$xmax
nad.bias.clust$ymin = 0
nad.bias.clust$ymax = sqrt(nad.bias.clust$integ/nad.bias.clust$integ.expected)
# NAD spie colors (cluster)
nad.clust.colors = rep(rev(brewer.pal(clust.cnt,"Blues")),2)
names(nad.clust.colors) = c('Act+','Act-','H3K9me','H3K27me','Lamin','Act+/NAD','Act-/NAD','H3K9me/NAD','H3K27me/NAD','Lamin/NAD')
nad.clust.lines = c(rev(brewer.pal(clust.cnt,"Blues")),rep('black',5))
names(nad.clust.lines) = c('Act+','Act-','H3K9me','H3K27me','Lamin','Act+/NAD','Act-/NAD','H3K9me/NAD','H3K27me/NAD','Lamin/NAD')
# Reorder cluster to avoid line overlap.
nad.bias.clust$clust = factor(nad.bias.clust$clust,c('Act+','Act-','H3K9me','H3K27me','Lamin','Act+/NAD','Act-/NAD','H3K9me/NAD','H3K27me/NAD','Lamin/NAD'))
nad.bias.clust = nad.bias.clust[order(nad.bias.clust$clust),]
# Generate spie with polar coordinates.
pdf(paste(fig.dir,'nup_clust_integ_spie.pdf',sep=''),useDingbats=F,width=7,height=7)
ggplot(nad.bias.clust) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax,color=clust),alpha=0.5) + coord_polar('x') + scale_fill_manual(name='3D pattern',values=nad.clust.colors) + scale_color_manual(name='3D pattern',values=nad.clust.lines) + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Integration enrichment in NAD (Nup153 Associated Domain) / 3D clusters')
dev.off()



# Cool pie chart.
## Generate minimal dataset.
clust.full$span.mb = round(clust.full$span/1e4)/1e2
cmini = clust.full[c('clust','cat','span.mb','integ.mb')]
## Reorder cluster/cat levels.
cmini$clust = factor(cmini$clust,c("Act+","Act-","H3K27me","H3K9me","Lamin"))
cmini$cat = as.factor(cmini$cat)
# Sort by cluster.
cmini = cmini[order(cmini$clust),]
## Compute y limits from cumsum of sizes.
cmini$ymax = cumsum(cmini$span.mb)
cmini$ymin = 0
cmini[2:nrow(cmini),]$ymin = cmini[1:(nrow(cmini)-1),]$ymax
# Plot, inner is cluster.
pdf(paste(fig.dir,'cluster_global_donut_cluster.pdf',sep=''),useDingbats=F,width=10,height=10)
ggplot(cmini) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=2.1,xmax=3)) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=3.1,xmax=4)) + scale_fill_manual(values = c(categ.colors,clust.colors)) + xlim(c(0,4)) + coord_polar('y') + theme_minimal() + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
dev.off()
# Sort by category.
cmini = cmini[order(cmini$cat),]
## Compute y limits from cumsum of sizes.
cmini$ymax = cumsum(cmini$span.mb)
cmini$ymin = 0
cmini[2:nrow(cmini),]$ymin = cmini[1:(nrow(cmini)-1),]$ymax
# Plot, inner is cluster.
pdf(paste(fig.dir,'cluster_global_donut_cat.pdf',sep=''),useDingbats=F,width=10,height=10)
ggplot(cmini) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=2.1,xmax=3)) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=3.1,xmax=4)) + scale_fill_manual(values = c(categ.colors,clust.colors)) + xlim(c(0,4)) + coord_polar('y') + theme_minimal() + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
dev.off()

# Plot integration density per compartment category.
pdf(paste(fig.dir,'global_cluster_hiv_density_perMb.pdf',sep=''),useDingbats=F,width=8,height=5)
ggplot(data=cmini,aes(fill=clust,x=cat,y=integ.mb))+geom_bar(stat='identity',position='dodge',colour='black') + xlab('Category') + ylab('HIV integrations per Mb') + scale_fill_manual(name='3D pattern',values=c(clust.colors))
dev.off()

# Expression boxplots
pdf(paste(fig.dir,'global_cluster_category_expression.pdf',sep=''),useDingbats=F,width=20,height=15)
#categ.expr$clust = as.factor(categ.expr$clust)
categ.expr$clust = factor(categ.expr$clust,c("Act+","Act-","H3K27me","H3K9me","Lamin"))
# Expression distribution by clusters.
p1 <- ggplot(categ.expr, aes(x = cat, y = expr, fill = cat)) + geom_boxplot(outlier.size=0.3) + scale_fill_manual(name='',values = categ.colors[clust.df$cat]) + coord_cartesian(ylim=c(0,3)) + theme_minimal() + xlab('') + ylab('Normalized HIV expression (log10)')
# Expression distribution by categories.
p2 <- ggplot(categ.expr, aes(x = clust, y = expr, fill = clust)) + geom_boxplot(outlier.size=0.3) + scale_fill_manual(name='',values = rev(brewer.pal(clust.cnt,"Blues"))) + coord_cartesian(ylim=c(0,3)) + theme_minimal() + xlab('') + ylab('')
# Expression by cluster and category (Grouped by cluster)
p3 <- ggplot(categ.expr, aes(x = clust, y = expr, fill = cat)) + geom_boxplot(outlier.size=0.3) + scale_fill_manual(name='Category',values = categ.colors[clust.df$cat]) + coord_cartesian(ylim=c(0,3)) + theme_minimal() + xlab('3D pattern') + ylab('Normalized HIV expression (log10)')
# Expression by cluster and category (Grouped by category)
p4 <- ggplot(categ.expr, aes(x = cat, y = expr, fill = clust)) + geom_boxplot(outlier.size=0.3) + scale_fill_manual(name='3D pattern',values = rev(brewer.pal(clust.cnt,"Blues"))) + coord_cartesian(ylim=c(0,3)) + theme_minimal() + xlab('Category') + ylab('')
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()

# Same figure for presentation
# Expression by cluster and category (Grouped by category)
pdf(paste(fig.dir,'global_cluster_category_expression_presentation.pdf',sep=''),useDingbats=F,width=4,height=4)
ggplot(categ.expr[categ.expr$cat == "actgen" | categ.expr$cat == "int",], aes(x = cat, y = expr, fill = clust)) + xlab("genome category") + ylab("HIV normalized expression (log)") + geom_boxplot(outlier.shape=NA) + scale_fill_manual(values = rev(brewer.pal(clust.cnt,"Blues"))) + coord_cartesian(ylim=c(0,2)) + theme_minimal()
dev.off()

# HIV Hotspots and clusters
# Make empty genome bins
genome$bins = ceiling(genome$size/hp.binsize)
hotspots = data.frame(chrom=c(),bin=c(),count=c())
for (i in seq(1,nrow(genome))) {
    hotspots = rbind(hotspots,data.frame(chrom=genome[i,1], bin=seq(0,as.integer(genome[i,3])-1), count=0))
}
# Read HIV integ site list
hp.all = read.table(hotspot.data)
names(hp.all) = c('chrom','locus')
hp.all$bin = floor(hp.all$locus/hp.binsize)
# Correct 2x read count for chrX and chrY
hp.all$count = 1
hp.all[hp.all$chrom=='chrX' | hp.all$chrom=='chrY',]$count = 2
hp.all = hp.all[,c('chrom','bin','count')]
# Append and aggregate counts
hotspots = rbind(hotspots,hp.all)
hotspots = aggregate(count ~ chrom+bin, data=hotspots, FUN=sum)
g.hp     = GRanges(Rle(hotspots$chrom),IRanges(start=hotspots$bin*hp.binsize,end=(hotspots$bin+1)*hp.binsize-1),
                         count=hotspots$count,cluster='NA')

# Compute 3D structure occupancy and expression of Top x%
hp.clust = data.frame(pctl=c(),cluster=c(),bin.count=c(),pctl.ratio=c())
exp.clust = data.frame(pctl=c(),expr=c())
hp.pctl  = c(0.99,0.95,0.9,(sum(g.hp$count==0)+1)/sum(g.hp$count>=0),0.0)
hp.names = c('Top 1%','Top 5%','Top 10%','All HIV','Genome')
for (i in seq(1,length(hp.pctl))) {
    p = hp.pctl[i]
    pctl.v = quantile(g.hp$count,p)
    hp.ov = findOverlaps(g.hp[g.hp$count >= pctl.v],clust.gr)
    hp.df = data.frame(pctl=hp.names[i],cluster=clust.gr[hp.ov@to]$global.cluster,count=1)
    hp.df = aggregate(count ~ pctl+cluster,data=hp.df,sum)
    hp.df$pctl.ratio = hp.df$count/sum(hp.df$count)
    hp.clust = rbind(hp.clust,hp.df)
    # BHIVE (expression) overlap
    exp.ov = findOverlaps(ghiv,g.hp[g.hp$count >= pctl.v])
    exp.clust = rbind(exp.clust,data.frame(pctl=hp.names[i],expr=ghiv[exp.ov@from]$expr))
}

# Reorder clusters.
hp.clust$cluster = factor(hp.clust$cluster,rev(c("Act+","Act-","H3K9me","H3K27me","Lamin")))
# Structure type: Stacked bar plot
pdf(paste(fig.dir,'hotspot_cluster_type_barplot.pdf',sep=''),useDingbats=F,width=6,height=8)
ggplot(hp.clust,aes(x=pctl,y=pctl.ratio,fill=cluster))+geom_bar(stat='identity') + scale_fill_manual(name='',values = brewer.pal(clust.cnt,"Blues")) + xlab('Hotspots') + ylab('3D structure type') + theme_minimal()
dev.off()

# HIV Expression: boxplots
# Remove 'Genome' category first:
exp.clust = exp.clust[exp.clust$pctl != 'Genome',]
pdf(paste(fig.dir,'hotspot_expression.pdf',sep=''),useDingbats=F,width=6,height=8)
ggplot(exp.clust,aes(x=pctl,y=expr,fill=pctl))+ geom_boxplot(outlier.shape=NA,alpha=0.3) + xlab('') + ylab('HIV normalized expression (log)') + ylim(c(0,2)) + theme_minimal() + guides(fill=guide_legend(title='Hotspots'))
dev.off()

# Compute 3D structure occupancy and expression of step hotspots 0-1%,1-5%,5-10%...
hp.dclust = data.frame(pctl=c(),cluster=c(),bin.count=c(),pctl.ratio=c())
exp.dclust = data.frame(pctl=c(),expr=c())
hp.topmax  = c(1.00,0.995,0.99,0.95,0.90,0.85,0.80)
hp.topmin  = c(0.995,0.99,0.95,0.90,0.85,0.80,sum(g.hp$count==0)/sum(g.hp$count>=0))
hp.dnames = c('Top [0,0.5)%','Top [0.5,1)%','Top [1,5)%','Top [5,10)%','Top [10,15)%','Top [15,20)%',paste('Top [20,',toString(as.integer(100-hp.topmin[7]*100)),')%',sep=''))
for (i in seq(1,length(hp.topmax))) {
    pctl.max = quantile(g.hp$count,hp.topmax[i])
    pctl.min = quantile(g.hp$count,hp.topmin[i])
    # Cluster overlap
    hp.ov = findOverlaps(g.hp[g.hp$count <= pctl.max & g.hp$count > pctl.min],clust.gr)
    hp.df = data.frame(pctl=hp.dnames[i],cluster=clust.gr[hp.ov@to]$global.cluster,count=1)
    hp.df = aggregate(count ~ pctl+cluster,data=hp.df,sum)
    hp.df$pctl.ratio = hp.df$count/sum(hp.df$count)
    hp.dclust = rbind(hp.dclust,hp.df)
    # BHIVE (expression) overlap
    exp.ov = findOverlaps(ghiv,g.hp[g.hp$count <= pctl.max & g.hp$count > pctl.min])
    exp.dclust = rbind(exp.dclust,data.frame(pctl=hp.dnames[i],expr=ghiv[exp.ov@from]$expr))
}

# Reorder clusters.
hp.dclust$cluster = factor(hp.dclust$cluster,rev(c("Act+","Act-","H3K9me","H3K27me","Lamin")))
# Structure type: Stacked bar plot
# Add sample size to labels
pdf(paste(fig.dir,'hotspot_cluster_type_step_barplot.pdf',sep=''),useDingbats=F,width=10,height=8)
ggplot(hp.dclust,aes(x=pctl,y=pctl.ratio,fill=cluster))+geom_bar(stat='identity') + scale_fill_manual(name='',values = brewer.pal(clust.cnt,"Blues")) + xlab('Genomic bins sorted by HIV integration density') + ylab('3D structure type') + theme_minimal()
dev.off()
# Hotspot expression: Box plots
exp.xlabels = c()
for (i in seq(1,length(hp.dnames))) {
    exp.xlabels = c(exp.xlabels,paste(hp.dnames[i],"\n\nn=",sum(exp.dclust$pctl==hp.dnames[i]),sep=''))
}
pdf(paste(fig.dir,'hotspot_step_expression.pdf',sep=''),useDingbats=F,width=8,height=8)
ggplot(exp.dclust,aes(x=pctl,y=expr,fill=pctl))+ geom_boxplot(outlier.shape=NA,alpha=0.3) + xlab('') + ylab('HIV normalized expression (log)') + ylim(c(0,2)) + theme_minimal() + guides(fill=guide_legend(title='Hotspots')) + scale_x_discrete(labels=exp.xlabels)
dev.off()

# Combined Hotspot 3D occupancy (Ribbon)
hp.cm = data.frame(pctl=c(),cluster=c(),bin.count=c(),pctl.ratio=c(),ymin=c(),ymax=c())
for (p in c(0.001,seq(0.01,0.45,0.005))) {
    pctl.v = quantile(g.hp$count,1-p)
    hp.ov = findOverlaps(g.hp[g.hp$count >= pctl.v],clust.gr)
    hp.df = data.frame(pctl=p*100,cluster=clust.gr[hp.ov@to]$global.cluster,count=1)
    # Add rows with count 0 to ensure all clusters are present for each x value!
    hp.df = rbind(hp.df,data.frame(pctl=p*100,cluster=c("Act+","Act-","H3K9me","H3K27me","Lamin"),count=0))
    hp.df = aggregate(count ~ pctl+cluster,data=hp.df,sum)
    hp.df$pctl.ratio = hp.df$count/sum(hp.df$count)
    hp.df$cluster = factor(hp.df$cluster,c("Act+","Act-","H3K9me","H3K27me","Lamin"))
    hp.df = hp.df[order(hp.df$cluster),]
    hp.df$ymax = cumsum(hp.df$pctl.ratio)
    hp.df$ymin = c(0,hp.df$ymax[1:nrow(hp.df)-1])
    hp.cm = rbind(hp.cm,hp.df)
}

pdf(paste(fig.dir,'hotspot_cluster_type_ribbon.pdf',sep=''),useDingbats=F,width=12,height=6)
ggplot(hp.cm) + geom_ribbon(aes(ymin=ymin,ymax=ymax,x=pctl,fill=cluster), alpha = 1) + xlim(c(0,40)) + theme_minimal() + xlab('Top x% bins with highest HIV density') + ylab('3D structure type') + scale_fill_manual(name='',values = rev(brewer.pal(clust.cnt,"Blues")))
dev.off()

# Podries fer un Ribbon tambe per expression!!
# A l'eix X hi poses el cluster density.
# A l'eix Y l'expressio i fas un plot continu amb la mediana (linia) i els percentils 10-90 (area). La part complicada es definir els limits dels grups de mostres, es a dir, si agafes de 0 a 0.5%, de 0.5% a 1%, etc... Perque si vols que sembli continua has d'agafar poques mostres pero alhora has d'assegurar-te que les distribucions tinguin suficients mostres. Tambe podries fer una sliding window, de mes a menys i que les distribucions estiguin repetides, encara que el que mostres a la X es el centre de la finestra!
#
# SLIDING WINDOW EM SEMBLA MOLT BONA IDEA PER FER EL RIBBON EXPRESSIO VS HOTSPOT HIV DENSITY.
#
