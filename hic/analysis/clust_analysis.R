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
# Also: set up percentile ranges and their names in section 6.1. (hp.pctl.intv)


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

# Analysis setup
pattern_sizes_pie    = TRUE
pattern_category_pie = TRUE
intrachromosomal_chromatin_heatmaps = TRUE
interchromosomal_chromatin_heatmap  = TRUE
interchromosomal_hiv_heatmap        = TRUE


# Function handlers
setdiff_nw      <- function(...) { suppressWarnings(setdiff(...)) };
findOverlaps_nw <- function(...) { suppressWarnings(findOverlaps(...)) };

# Helper functions
MaxTable <- function(InVec) {
  if (!is.factor(InVec)) InVec <- factor(InVec)
  t       <- table(InVec)
  maxVals <- names(t)[t == max(t)]
  if (length(maxVals) > 1)
      NA
  else
      maxVals
}


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
## 2.2. 3D patterns.
## 2.3. Category background of 3D patterns
##
## 3. Signal distribution in 3D patterns.
## 3.1. Chromatin marks in 3D patterns.
## 3.2. HIV in 3D patterns.
##
## 4. HIV integration pattern
## 4.1. Integration density in Categories.
## 4.2. Integration density in 3D patterns.
## 4.3. Local integration density in categories (grouped and normalized by 3D pattern)
## 4.5. HIV integration density (grouped by 3D pattern or category)
## 
## 5. Gene expression in 3D patterns
## 5.1. Endogenous expression in 3D patterns
## 5.2. HIV expression in 3D patterns
##
## 6. HIV hotspots
## 6.1. Hotspot occupancy
## 6.2. HIV expression in hotspots.


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
genome      <- read.table(chrsz.data, as.is=TRUE, col.names=c('chr','size'))
genome.gr   <- GRanges(Rle(genome$chr), IRanges(start=1, width=genome$size))
genome$bins <- ceiling(genome$size/hp.binsize)
genome.bins <- data.frame(chrom=rep(genome$chr,genome$bins), bin=unlist(lapply(genome$bins, seq))-1, count=0)
# Verbose genome
cat('[1.4.1] chromosomes: ',length(genome.gr),'\n',sep='')
cat('[1.4.1] coverage: ',sum(as.numeric(sum(coverage(genome.gr)))),' nt\n',sep='')
cat('[1.4.1] hotspot bin size: ',hp.binsize,'\n',sep='')
cat('[1.4.1] hotspot bin count: ',sum(genome$bins),'\n',sep='')

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
## 2.2. 3D patterns.
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
## 2.3. Category background of 3D patterns
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


##
## 4.1. Integration density in Categories.
##
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


##
## 4.2. Integration density in 3D patterns.
##
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


##
## 4.3. Local integration density in categories (grouped and normalized by 3D pattern)
##
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


##
## 4.4. Local integration density in 3D patterns (grouped and normalized by Category)
##
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


##
## 4.5. HIV integration density (grouped by 3D pattern or category)
##
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


##
## 5.1 Endogenous expression in 3D patterns
##
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


##
## 5.2 HIV expression in 3D patterns
##
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

# Boxplot figures.
figure.fn <- paste(fig.dir,'expr_bhive.pdf',sep='')
cat('[5.2] HIV expression in 3D patterns (full figure) -> ',figure.fn,'\n',sep='')
pdf(figure.fn,useDingbats=F,width=20,height=15)
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()


##
## 6. HIV hotspots
##
cat('[6. HIV hotspots]\n')

# Read HIV integ site list
hp.all       <- read.table(hotspot.data, col.names= c('chrom','locus'))
hp.all$bin   <- floor(hp.all$locus/hp.binsize)
# Correct 2x read count for (chrX, chrY)
hp.all$count                                               <- 1
hp.all[hp.all$chrom=='chrX' | hp.all$chrom=='chrY',]$count <- 2
# Append and aggregate counts
hotspots      <- aggregate(count ~ chrom+bin, data=rbind(genome.bins, hp.all[,c('chrom','bin','count')]), FUN=sum)
# Compute hotspot percentiles.
hotspots$pctl <- ecdf(hotspots$count)(hotspots$count)

# Define percentile ranges.
hp.pctl.intv    <- c(0.995,0.99,0.95,0.9,(sum(hotspots$count==0)+1)/sum(hotspots$count>=0),0.0)
hp.pctl.names   <- c('Top 0.5%','Top 0.5%-1%','Top 1%-5%','Top 5%-10%','Sparse HIV','No HIV')
# Sort intervals
hp.pctl.ord     <- order(hp.pctl.intv)
hp.pctl.intv.s  <- hp.pctl.intv[hp.pctl.ord]
hp.pctl.names.s <- hp.pctl.names[hp.pctl.ord]

# Make GenomicRanges and intersect with cluster/category information.
hp.gr         <- GRanges(seqnames = Rle(hotspots$chrom),
                     ranges   = IRanges(start=hotspots$bin*hp.binsize,end=(hotspots$bin+1)*hp.binsize-1),
                     count    = hotspots$count,
                     pctl     = hotspots$pctl,
                     cluster  = NA,
                     intv     = NA
                     )

##
## 6.1. Hotspot occupancy
##
cat('[6.1. Hotspot occupancy in 3D patterns]\n')

# Compute predominant cluster in hotspot bins
hp.clu.ov     <- findOverlaps(hp.gr, clust.categ.gr, ignore.strand = TRUE)
hp.clu.ref    <- data.frame(hp.bin = hp.clu.ov@from, clu = clust.categ.gr[hp.clu.ov@to]$clust)
hp.clu.ref    <- aggregate(clu ~ hp.bin, data = hp.clu.ref, FUN=MaxTable)

hp.gr$cluster[hp.clu.ref$hp.bin] <- hp.clu.ref$clu

# Classify bins in percentile intervals.
hp.gr$intv    <- factor(hp.pctl.names.s[findInterval(hp.gr$pctl,hp.pctl.intv.s)], hp.pctl.names)

# Compute 3D pattern rate for all intervals.
hp.clu.table  <- table(
    data.frame(
        intv  = hp.gr$intv,
        clust = factor(hp.gr$cluster,rev(c("Act+","Act-","H3K27me","H3K9me","Lamin")))
    )
)
hp.clu.rates  <- as.data.frame(hp.clu.table/rowSums(hp.clu.table))

# Structure type: Stacked bar plot
figure.fn <- paste(fig.dir,'hiv_hotspot_3dpatterns.pdf',sep='')
cat('[6.1] Distribution of HIV hotspots in 3D patterns -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(hp.clu.rates, aes(x=intv, y=Freq, fill=clust)) +
    geom_bar(stat='identity') +
    scale_fill_manual(name='', values = clust.colors) +
    xlab(paste('Genomic bins (',round(hp.binsize/1e3),' kb) sorted by HIV density',sep='')) +
    ylab('3D pattern ratio') +
    theme_minimal()
dev.off()


##
## 6.2. HIV expression in hotspots.
##
cat('[6.2. HIV expression in hotspots]\n')

# Cross expression data with hotspot/cluster information.
hp.clu.expr.ov <- findOverlaps_nw(bhive.gr, hp.gr, ignore.strand = TRUE)
hp.clu.expr    <- data.frame(
    intv  = hp.gr$intv[hp.clu.expr.ov@to],
    clust = factor(hp.gr$cluster[hp.clu.expr.ov@to], c("Act+","Act-","H3K27me","H3K9me","Lamin")),
    expr  = bhive.gr$expr[hp.clu.expr.ov@from]
)

# Expression of hotspot intervals
figure.fn <- paste(fig.dir,'hiv_expr_hotspot.pdf',sep='')
cat('[6.2] Expression in hotspots -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=8, height=8)
ggplot(hp.clu.expr,aes(x=intv,y=expr,fill=intv)) +
    geom_boxplot(outlier.shape=NA, alpha=0.3) +
    xlab('') +
    ylab('HIV normalized expression (log)') +
    ylim(c(0,2)) +
    theme_minimal() +
    guides(fill=guide_legend(title='Hotspots'))
dev.off()

# Expression of hotspot intervals, grouped by 3d pattern (Act+, Act- and Lamin only)
figure.fn <- paste(fig.dir,'hiv_expr_hotspot_3dpattern.pdf',sep='')
cat('[6.2] Expression in hotspots, grouped by 3D Pattern -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=12, height=8)
ggplot(hp.clu.expr[hp.clu.expr$clust %in% c('Act+','Act-','Lamin'),],aes(x=clust,y=expr,fill=intv)) +
    geom_boxplot(outlier.shape=NA, alpha=0.3) +
    xlab('') +
    ylab('HIV normalized expression (log)') +
    ylim(c(0,2)) +
    theme_minimal() +
    guides(fill=guide_legend(title='Hotspots'))
dev.off()


# Remove H3K27me (almost no HIV), also remove H3K9me from Top 0.5-1% (poor stats)
hp.clu.tmp <- hp.clu.expr[hp.clu.expr$clust %in% c('Act+','Act-','H3K9me','Lamin'),]
hp.clu.tmp <- hp.clu.tmp[!(hp.clu.tmp$clust=='H3K9me' & hp.clu.tmp$intv=='Top 0.5%-1%'),]

# Expression in 3d patterns grouped by interval (The more sparse, the greater the expression difference)
figure.fn <- paste(fig.dir,'hiv_expr_3dpattern_hotspot.pdf',sep='')
cat('[6.2] Expression in 3D patterns, grouped by hotspots -> ',figure.fn,'\n',sep='')
pdf(figure.fn, useDingbats=F, width=12, height=8)
ggplot(hp.clu.tmp, aes(x=intv, y=expr, fill=clust)) +
    geom_boxplot(outlier.shape=NA, alpha=1) +
    scale_fill_manual(values=clust.colors, name='3D Pattern') +
    xlab('') +
    ylab('HIV normalized expression (log)') +
    ylim(c(0,2)) +
    theme_minimal()
dev.off()


##
## End of pipeline
## 

##
## OTHER ANALYSES
##

## #   - Load Nup data.
#nup.damid    = 'chip_datasets/Nup153_IMR90.bed'
## nup.data = read.table(nup.damid)
## gnup = GRanges(seqnames=Rle(nup.data$V1),ranges=IRanges(start=nup.data$V2,end=nup.data$V3))

## # NAD categories
## nad.categ = list()
## nad.categ[['silgen-nad']] = suppressWarnings(intersect(categ[['silgen']],gnup))
## nad.categ[['actgen-nad']] = suppressWarnings(intersect(categ[['actgen']],gnup))
## nad.categ[['actpro-nad']] = suppressWarnings(intersect(categ[['actpro']],gnup))
## nad.categ[['enh-nad']] = suppressWarnings(intersect(categ[['enh']],gnup))
## nad.categ[['int-nad']] = suppressWarnings(intersect(categ[['int']],gnup))
## # Non-NAD categories
## nad.categ[['silgen']] = setdiff_nw(categ[['silgen']],gnup)
## nad.categ[['actgen']] = setdiff_nw(categ[['actgen']],gnup)
## nad.categ[['actpro']] = setdiff_nw(categ[['actpro']],gnup)
## nad.categ[['enh']] = setdiff_nw(categ[['enh']],gnup)
## nad.categ[['int']] = setdiff_nw(categ[['int']],gnup)

## # SuperSpyChart Nup153
## #   - Add category info to clusters.
## clust.gr$cat = ''
## #   * The order of the loop is important to correctly resolve bins falling in category edges.
## for (idx in c('int','silgen','actgen','actpro','enh')) {
##     ov = suppressWarnings(findOverlaps(clust.gr,categ[[idx]]))
##     clust.gr[ov@from]$cat = idx
## }

## nup.all = data.frame(clust=c(),cat=c(),bg.nt=c(),nup.nt=c())
## for (i in unique(clust.gr$global.cluster)) {
##     for (idx in names(categ)) {
##         gc = clust.gr[clust.gr$global.cluster == i & clust.gr$cat == idx]
##         nup.df = data.frame(clust=toString(i),cat=idx,bg.nt=sum(as.numeric(sum(coverage(gc)))),nup.nt=sum(as.numeric(sum(coverage(intersect(gc,gnup,ignore.strand=T))))))
##         nup.df$cat = levels(droplevels(nup.df$cat))
##         nup.all = rbind(nup.all,nup.df)
##     }
## }

## #   Spie on category.
## nup.cat = merge(aggregate(bg.nt ~ cat,data=nup.all,FUN=sum),aggregate(nup.nt ~ cat,data=nup.all,FUN=sum),by='cat')
## nup.cat$cat = factor(nup.cat$cat,c("actgen","actpro","enh","silgen","int"))
## nup.cat = nup.cat[order(nup.cat$cat),]
## nup.cat$xmax = cumsum(nup.cat$bg.nt)
## nup.cat$xmin = c(0,nup.cat[1:(nrow(nup.cat)-1),]$xmax)
## nup.cat$ymin = 0
## nup.cat$ymax = sqrt(nup.cat$nup.nt/sum(nup.cat$nup.nt)/(nup.cat$bg.nt/sum(nup.cat$bg.nt)))
## pdf(paste(fig.dir,'nup_category_signal_spie.pdf',sep=''),useDingbats=F,width=8,height=8)
## ggplot(nup.cat) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color=categ.colors[nup.cat$cat]) + scale_fill_manual(name='Category',values=categ.colors)+ coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Nup153 enrichment')
## dev.off()
## #   Spie on 3D pattern.
## nup.clu = merge(aggregate(bg.nt ~ clust,data=nup.all,FUN=sum),aggregate(nup.nt ~ clust,data=nup.all,FUN=sum),by='clust')
## nup.clu$clust = factor(nup.clu$clust,c("Act+","Act-","H3K27me","H3K9me","Lamin"))
## nup.clu = nup.clu[order(nup.clu$clust),]
## nup.clu$xmax = cumsum(nup.clu$bg.nt)
## nup.clu$xmin = 0
## nup.clu[2:nrow(nup.clu),]$xmin = nup.clu[1:(nrow(nup.clu)-1),]$xmax
## nup.clu$ymin = 0
## nup.clu$ymax = sqrt(nup.clu$nup.nt/sum(nup.clu$nup.nt)/(nup.clu$bg.nt/sum(nup.clu$bg.nt)))
## pdf(paste(fig.dir,'nup_clust_signal_spie.pdf',sep=''),useDingbats=F,width=8,height=8)
## ggplot(nup.clu) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5,color="#555555") + scale_fill_manual(name='3D pattern',values=clust.colors)+ coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Nup153 enrichment')
## dev.off()

## # Crazy barplot
## nup.all$clust = factor(nup.all$clust,c("Act+","Act-","H3K27me","H3K9me","Lamin"))
## nup.all$cat = factor(nup.all$cat,c("actgen","actpro","enh","silgen","int"))
## nup.all = nup.all[order(nup.all$cat,nup.all$clust),]
## nup.all$ratio.bg = nup.all$bg.nt / sum(nup.all$bg.nt)
## nup.all$ratio.nup = nup.all$nup.nt / sum(nup.all$nup.nt)
## nup.all$xmax = cumsum(nup.all$bg.nt)
## nup.all$xmin = 0
## nup.all[2:nrow(nup.all),]$xmin = nup.all[1:(nrow(nup.all)-1),]$xmax
## nup.all$ymax = sqrt(nup.all$ratio.nup/nup.all$ratio.bg)
## nup.all$ymin = 0

## nup.sq = nup.cat[,c('cat','xmin','xmax','ymin','ymax')]
## nup.sq$clust = nup.sq$cat
## for (c in unique(nup.sq$cat)) {
##     nup.all[nup.all$cat == c,'ymin'] = nup.sq[nup.sq$cat == c,'ymax']
## }
## nup.sq = rbind(nup.sq,nup.all[,c('cat','xmin','xmax','ymin','ymax','clust')])
## clust.colors.gray = rep("#555555",5)
## names(clust.colors.gray) = names(clust.colors)
## pdf(paste(fig.dir,'nup_clust_cat_signal_superspie.pdf',sep=''),useDingbats=F,width=12,height=12)
## ggplot(nup.sq) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax,color=clust),alpha=0.5) + scale_fill_manual(name='',values=c(clust.colors,categ.colors)) + scale_color_manual(name='',values=c(clust.colors.gray,categ.colors))+ coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Nup153 enrichment by Category and 3D pattern')
## dev.off()

## # Nup Barplot all. (it does not reflect well the reality beacuse the bars may represent tiny or huge regions and are all the same width).
## nup.all$cov = nup.all$nup.nt/nup.all$bg.nt
## ggplot(nup.all,aes(fill=clust,x=cat,y=cov)) + geom_bar(stat='identity',position="dodge",color="#555555") + scale_fill_manual(values=clust.colors)

## # NAD integration bias vs category.
## nad.bias.cat = data.frame(cat=c(),span=c(),cat.ratio=c(),integ=c(),integ.expected=c(),integ.mb=c())
## total.span = sum(clust.categ$span)
## for (i in names(nad.categ)) {
##     gcat = nad.categ[[i]]
##     ov.cnt = suppressWarnings(countOverlaps(gcat,ghiv,ignore.strand=T))
##     span = sum(as.numeric(sum(coverage(gcat))))
##     intexpect = nrow(hiv) * span/total.span
##     integ = sum(ov.cnt)
##     nad.bias.cat = rbind(nad.bias.cat, data.frame(cat=i,span=span,cat.ratio=span/total.span,integ=integ,integ.expected=intexpect, integ.mb=integ*1e6/span))
## }
## # Prepare data (xmin,xmax,ymin,ymax) for spie chart.
## nad.bias.cat$cat = factor(nad.bias.cat$cat,c("actgen/NAD","actgen","actpro/NAD","actpro","enh/NAD","enh","silgen/NAD","silgen","int/NAD","int"))
## nad.bias.cat = nad.bias.cat[order(nad.bias.cat$cat),]
## nad.bias.cat$xmax = cumsum(nad.bias.cat$span)
## nad.bias.cat$xmin = 0
## nad.bias.cat[2:nrow(nad.bias.cat),]$xmin = nad.bias.cat[1:(nrow(nad.bias.cat)-1),]$xmax
## nad.bias.cat$ymin = 0
## nad.bias.cat$ymax = sqrt(nad.bias.cat$integ/nad.bias.cat$integ.expected)

## # NAD spie colors (category)
## nad.cat.colors = rep(c("#92C46DA0","#548B54","#000000A0", "#19485780", "#C0955FE0"),2)
## names(nad.cat.colors) = c('actgen','actpro','silgen','enh','int','actgen/NAD','actpro/NAD','silgen/NAD','enh/NAD','int/NAD')
## nad.cat.lines = c(rep('blank',5),rep('solid',5))
## names(nad.cat.lines) = c('actgen','actpro','silgen','enh','int','actgen/NAD','actpro/NAD','silgen/NAD','enh/NAD','int/NAD')

## # Generate spie with polar coordinates.
## pdf(paste(fig.dir,'nup_category_integ_spie.pdf',sep=''),useDingbats=F,width=7,height=7)
## ggplot(nad.bias.cat) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),alpha=0.5) + coord_polar('x') + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Integration enrichment in NAD (Nup153 Associated Domain) / Categories') + scale_linetype_manual(name='Category',values=nad.cat.lines) + scale_fill_manual(name='Category',values=nad.cat.colors)
## dev.off()

## # NAD integration bias vs 3d pattern.
## nad.bias.clust = data.frame(clust=c(),span=c(),clust.ratio=c(),integ=c(),integ.expected=c(),integ.mb=c())
## for (i in unique(clust.gr$global.cluster)) {
##     gclu = clust.gr[clust.gr$global.cluster == toString(i),]
##     gclu.nad = suppressWarnings(intersect(gclu,gnup,ignore.strand=T))
##     gclu.non = setdiff_nw(gclu,gnup,ignore.strand=T)
##     nad.cnt = sum(suppressWarnings(countOverlaps(gclu.nad,ghiv,ignore.strand=T)))
##     non.cnt = sum(suppressWarnings(countOverlaps(gclu.non,ghiv,ignore.strand=T)))
##     nad.span = sum(as.numeric(sum(coverage(gclu.nad))))
##     non.span = sum(as.numeric(sum(coverage(gclu.non))))
##     nad.intexpect = nrow(hiv) * nad.span/total.span
##     non.intexpect = nrow(hiv) * non.span/total.span
##     nad.bias.clust = rbind(nad.bias.clust, data.frame(clust=paste(i,"/NAD",sep=""),span=nad.span,clust.ratio=nad.span/total.span,integ=nad.cnt,integ.expected=nad.intexpect, integ.mb=nad.cnt*1e6/span))
##     nad.bias.clust = rbind(nad.bias.clust, data.frame(clust=i,span=non.span,clust.ratio=non.span/total.span,integ=non.cnt,integ.expected=non.intexpect, integ.mb=non.cnt*1e6/span))
## }
## # Prepare data (xmin,xmax,ymin,ymax) for spie chart.
## nad.bias.clust$clust = factor(nad.bias.clust$clust,c("Act+/NAD","Act+","Act-/NAD","Act-","H3K27me/NAD","H3K27me","H3K9me/NAD","H3K9me","Lamin/NAD","Lamin"))
## nad.bias.clust = nad.bias.clust[order(nad.bias.clust$clust),]
## nad.bias.clust$xmax = cumsum(nad.bias.clust$span)
## nad.bias.clust$xmin = 0
## nad.bias.clust[2:nrow(nad.bias.clust),]$xmin = nad.bias.clust[1:(nrow(nad.bias.clust)-1),]$xmax
## nad.bias.clust$ymin = 0
## nad.bias.clust$ymax = sqrt(nad.bias.clust$integ/nad.bias.clust$integ.expected)
## # NAD spie colors (cluster)
## nad.clust.colors = rep(rev(brewer.pal(clust.cnt,"Blues")),2)
## names(nad.clust.colors) = c('Act+','Act-','H3K9me','H3K27me','Lamin','Act+/NAD','Act-/NAD','H3K9me/NAD','H3K27me/NAD','Lamin/NAD')
## nad.clust.lines = c(rev(brewer.pal(clust.cnt,"Blues")),rep('black',5))
## names(nad.clust.lines) = c('Act+','Act-','H3K9me','H3K27me','Lamin','Act+/NAD','Act-/NAD','H3K9me/NAD','H3K27me/NAD','Lamin/NAD')
## # Reorder cluster to avoid line overlap.
## nad.bias.clust$clust = factor(nad.bias.clust$clust,c('Act+','Act-','H3K9me','H3K27me','Lamin','Act+/NAD','Act-/NAD','H3K9me/NAD','H3K27me/NAD','Lamin/NAD'))
## nad.bias.clust = nad.bias.clust[order(nad.bias.clust$clust),]
## # Generate spie with polar coordinates.
## pdf(paste(fig.dir,'nup_clust_integ_spie.pdf',sep=''),useDingbats=F,width=7,height=7)
## ggplot(nad.bias.clust) + geom_hline(yintercept=sqrt(c(0.5,2,3)),alpha=0.5,linetype="dotted") + geom_hline(yintercept=sqrt(c(1)),alpha=0.5) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax,color=clust),alpha=0.5) + coord_polar('x') + scale_fill_manual(name='3D pattern',values=nad.clust.colors) + scale_color_manual(name='3D pattern',values=nad.clust.lines) + theme_minimal() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + ggtitle('Integration enrichment in NAD (Nup153 Associated Domain) / 3D clusters')
## dev.off()

## # Cool pie chart.
## ## Generate minimal dataset.
## clust.full$span.mb = round(clust.full$span/1e4)/1e2
## cmini = clust.full[c('clust','cat','span.mb','integ.mb')]
## ## Reorder cluster/cat levels.
## cmini$clust = factor(cmini$clust,c("Act+","Act-","H3K27me","H3K9me","Lamin"))
## cmini$cat = as.factor(cmini$cat)
## # Sort by cluster.
## cmini = cmini[order(cmini$clust),]
## ## Compute y limits from cumsum of sizes.
## cmini$ymax = cumsum(cmini$span.mb)
## cmini$ymin = 0
## cmini[2:nrow(cmini),]$ymin = cmini[1:(nrow(cmini)-1),]$ymax
## # Plot, inner is cluster.
## pdf(paste(fig.dir,'cluster_global_donut_cluster.pdf',sep=''),useDingbats=F,width=10,height=10)
## ggplot(cmini) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=2.1,xmax=3)) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=3.1,xmax=4)) + scale_fill_manual(values = c(categ.colors,clust.colors)) + xlim(c(0,4)) + coord_polar('y') + theme_minimal() + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
## dev.off()
## # Sort by category.
## cmini = cmini[order(cmini$cat),]
## ## Compute y limits from cumsum of sizes.
## cmini$ymax = cumsum(cmini$span.mb)
## cmini$ymin = 0
## cmini[2:nrow(cmini),]$ymin = cmini[1:(nrow(cmini)-1),]$ymax
## # Plot, inner is cluster.
## pdf(paste(fig.dir,'cluster_global_donut_cat.pdf',sep=''),useDingbats=F,width=10,height=10)
## ggplot(cmini) + geom_rect(aes(fill=cat,ymin=ymin,ymax=ymax,xmin=2.1,xmax=3)) + geom_rect(aes(fill=clust,ymin=ymin,ymax=ymax,xmin=3.1,xmax=4)) + scale_fill_manual(values = c(categ.colors,clust.colors)) + xlim(c(0,4)) + coord_polar('y') + theme_minimal() + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
## dev.off()

