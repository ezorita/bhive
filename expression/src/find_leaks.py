import sys

f    = open(sys.argv[1])
last = ''
reps = list()

for line in f:
    l    = line.split('\t')
    loc = l[2]+l[3]+l[4]
    if loc == last:
        reps.append(line)
    else:
        if len(reps) > 1:
            for b in reps:
                print b.rstrip()
        del reps[:]
        reps.append(line)
    last = loc
        
if len(reps) > 1:
    for b in reps:
        print b
