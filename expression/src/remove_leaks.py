import sys

f    = open(sys.argv[1])
last = ''
reps = list()

# Read header
header = f.readline()
print header.rstrip()
header = header.rstrip().split('\t')
reads_idx = header.index('reads')
dna_idx = header.index('dna')
rna_idx = header.index('rna')

for line in f:
    l    = line.split('\t')
    loc = l[2]+l[3]+l[4]
    if loc == last:
        reps.append(line)
    else:
        if len(reps) > 1:
            score = 0
            for b in reps:
                # [7] DNA count is NA
                # [8] RNA count is NA
                l = b.split('\t')
                if not (l[dna_idx] == 'NA' or l[dna_idx] == 0 or l[rna_idx] == 'NA'):
                   print b.rstrip()
        else:
            if len(reps) > 0:
                print reps[0].rstrip()
                
        del reps[:]
        reps.append(line)
    last = loc
        
if len(reps) > 0:
    if len(reps) > 1:
        score = 0
        for b in reps:
           # [7] DNA count is NA
           # [8] RNA count is NA
            l = b.split('\t')
            if not (l[dna_idx] == 'NA' or l[dna_idx] == 0 or l[rna_idx] == 'NA'):
                print b.rstrip()
    else:
        print reps[0].rstrip()
