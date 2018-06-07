import sys

f    = open(sys.argv[1])
last = ''
reps = list()

for line in f:
    l    = line.split('\t')
    brcd = l[0]
    if brcd == last:
        reps.append(line)
    else:
        if len(reps) > 1:
            for b in reps:
                print b.rstrip()
        del reps[:]
        reps.append(line)
    last = brcd
        
if len(reps) > 1:
    for b in reps:
        print b
