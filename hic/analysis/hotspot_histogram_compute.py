import numpy as np
import cooler

bin_factor = 20 #100kb bins
hic_bin    = 5e3
bin_size   = int(hic_bin*bin_factor)

percentiles = [[99,1],[95,10],[90,10],[0,10]]
max_samp_size = 200
sample_seed = 1
contrl_seed = 2
integ_file = 'integ_files/jurkat_bhive_bushman_chronly_sorted.integ'
chrsz_file = 'annotations/hg19_chromsizes.txt'

# Prealloc genome bins
integs = {}
with open(chrsz_file,'r') as f:
    for line in f:
        l = line.rstrip().split('\t')
        integs[l[0]] = np.zeros(np.int(l[1])/bin_size+1)

# Load integration files
with open(integ_file,'r') as f:
    for line in f:
        l = line.rstrip().split('\t')
        chrom = l[0]
        locus = int(l[1])
        if not integs.has_key(chrom):
            continue
        if chrom == 'chrX' or chrom == 'chrY':
            integs[chrom][locus/bin_size] += 2
        else:
            integs[chrom][locus/bin_size] += 1

# Make a big bin count list
allinteg = np.array([])
for chrom,arr in integs.iteritems():
    allinteg = np.append(allinteg,arr)

# Initialize random seed
np.random.seed(sample_seed)

## HOTSPOT SAMPLE
# Compute hotspot percentile
samp_size = {}
max_pctl = np.max([x[0] for x in percentiles])
max_v = np.percentile(allinteg,max_pctl)
# Compute sample size from highest percentile
for chrom in integs:
    samp_size[chrom] = min(len(np.where(integs[chrom] >= max_v)[0]),max_samp_size)
# Make hotspots samples
hotspots = {}
hot_samp = {}
contacts = {}
for [pctl,n_samp] in percentiles:
    pctl_v = np.percentile(allinteg,pctl)
    hotspots[pctl] = {}
    hot_samp[pctl] = {}
    contacts[pctl] = [[] for x in xrange(0,n_samp)]
    for chrom in integs:
        hotspots[pctl][chrom] = np.where(integs[chrom] >= pctl_v)[0]
        if samp_size[chrom] == 0: continue
        for s in xrange(0,n_samp):
            samp = np.random.choice(hotspots[pctl][chrom],size=samp_size[chrom],replace=False)
            if not hot_samp[pctl].has_key(chrom):
                hot_samp[pctl][chrom] = [samp,]
            else:
                hot_samp[pctl][chrom].append(samp)


# Load Cooler file.
cf = cooler.Cooler('cool_files/WT_hg19_5k_q10.cool')
hs0 = hot_samp[max_pctl]
# Compute all-pairs signal
for a in xrange(0,len(hs0)):
    chr_a  = hs0.keys()[a]
    for b in xrange(a+1,len(hs0)):
        chr_b  = hs0.keys()[b]
        print "{},{}".format(chr_a,chr_b)
        # Load Hi-C interchromosomal matrix
        m = cf.matrix(balance=False).fetch(chr_a,chr_b)
        # Compute pair-wise contact signal for hotspots.
        for [pctl,samp_cnt] in percentiles:
            for samp_no in xrange(0,samp_cnt):
                for i in hot_samp[pctl][chr_a][samp_no]:
                    for j in hot_samp[pctl][chr_b][samp_no]:
                        contacts[pctl][samp_no].append(np.sum(m[i*bin_factor:(i+1)*bin_factor,:][:,j*bin_factor:(j+1)*bin_factor]))

# Write output file
with open('hotspot_contact_histograms.txt','w+') as fout:
    for [pctl,n_samp] in percentiles:
        for samp_no in xrange(0,n_samp):
            hcont = np.array(contacts[pctl][samp_no])
            [hcont_hist,hcont_break] = np.histogram(hcont,bins=np.max(hcont))
            for (x,y) in zip(hcont_break,hcont_hist):
                fout.write("{}\t{}\t{}\t{}\n".format(pctl,samp_no,x,y))
