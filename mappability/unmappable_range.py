import sys

if len(sys.argv) != 3:
    sys.stdout.write("usage: "+sys.argv[0]+" <min_Q> <mapped_genome.sam>\n")
    sys.exit(1)

Q = int(sys.argv[1])
file = sys.argv[2]
cur_chr = ""
cur_beg = 0
cur_end = 0
cur_maq = 1

with open(file) as f:
    for line in f:
        # Ignore comments
        if line[0] == '#': continue
        # Split line
        map = line.rstrip().split('\t')
        # Parse SAM line
        if len(map) < 5: continue
        locus = map[0].split('_')
        chr = '_'.join(locus[:-2])
        beg = int(locus[-2])
        end = int(locus[-1])
        maq = int(map[4]) >= Q
        # Chromosome boundary
        if cur_chr != chr:
            # Print last interval
            if cur_maq == False:
                sys.stdout.write(chr+"\t"+str(cur_beg)+"\t"+str(cur_end)+"\n")
            cur_maq = 1
            cur_chr = chr
            
        # Intrachromosome interval
        if cur_maq == True and maq == False:
            cur_beg = beg
            cur_end = end
        elif cur_maq == False and maq == True:
            sys.stdout.write(chr+"\t"+str(cur_beg)+"\t"+str(beg-1)+"\n")
        elif cur_maq == False and maq == False:
            cur_end = end
        # Update region type
        cur_maq = maq

# Last interval before EOF
if cur_maq == False:
    sys.stdout.write(chr+"\t"+str(cur_beg)+"\t"+str(cur_end)+"\n")
