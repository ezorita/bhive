import sys

if len(sys.argv) != 4:
    sys.stderr.write("usage: "+sys.argv[0]+" <kmer size> <offset> <genome.fasta>\n")
    sys.exit(1)

k = int(sys.argv[1])
o = int(sys.argv[2])
file = sys.argv[3]
chr = "unknown_chromosome"
buffer = ""

with open(file) as f:
    for line in f:
        if line[0] == '>':
            for i in xrange(0,(len(buffer)-k)/o+1):
                sys.stdout.write(">"+chr+"_"+str(o*i+1)+"_"+str(o*i+k)+"\n"+buffer[o*i:(o*i+k)]+"\n")
            if len(buffer) >= k and (len(buffer)-k)%o > 0:
                sys.stdout.write(">"+chr+"_"+str(len(buffer)-k+1)+"_"+str(len(buffer))+"\n"+buffer[len(buffer)-k:]+"\n")
            chr = line[1:].rstrip()
            buffer = ""
        else:
            buffer += line.rstrip()

# Final print
for i in xrange(0,(len(buffer)-k)/o+1):
    sys.stdout.write(">"+chr+"_"+str(o*i+1)+"_"+str(o*i+k)+"\n"+buffer[o*i:(o*i+k)]+"\n")
if len(buffer) >= k and (len(buffer)-k)%o > 0:
    sys.stdout.write(">"+chr+"_"+str(len(buffer)-k+1)+"_"+str(len(buffer))+"\n"+buffer[len(buffer)-k:]+"\n")
