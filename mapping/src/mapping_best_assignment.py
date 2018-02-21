import argparse
import sys
import math
import re

parser = argparse.ArgumentParser(description='Find best barcode-locus assignment.')
parser.add_argument('-q', required=False, default=20, type=int, help="minimum mapping score (default 20)")
parser.add_argument('--min-reads', required=False,type=int,default=10,help="minimum # of reads to validate an assignment (default 10)")
parser.add_argument('--min-score', required=False,type=int,default=10,help="minimum assignment score (default 10)")
parser.add_argument('mapping_sam', help="SAM alignment file of integration loci with barcode in sequence name")
parser.add_argument('integ_output',help="Output file for valid integrations.")
parser.add_argument('filter_output',help="Output file for filtered/recombinant integrations.")
params = parser.parse_args()

bwa_q_thr = params.q
min_reads = params.min_reads
min_score = params.min_score
loc_merge_dist = 100 # Loci with same barcode at distance <= loc_merge_dist [bp] will be considered the same insertion
strand = ['-','+']
maps = {}
bcdset = set()

integ_file = open(params.integ_output,"w")
recomb_file = open(params.filter_output,"w")

with open(params.mapping_sam) as f:
   for line in f:
      if line[0] == '@': continue
      m = line.split('\t')
      # Exclude self-ligated reads.
      if m[2] == '*' or int(m[4]) < bwa_q_thr: continue
      # Parse SAM fields
      barcode = m[0]
      rchr = m[2]
      rloc = int(m[3])
      flags = int(m[1])
      dirflag = (flags & 16) > 0
      cigar = m[5]

      # Don't allow soft clipping of > 3nt at the LTR end.
      asymb = re.split('[0-9]+',cigar)[1:]
      avals = re.split('M|S|I|H|D',cigar)[:-1]
      if (dirflag and (asymb[-1] in ['S','H']) and int(avals[-1]) > 3) or (not dirflag and (asymb[0] in ['S','H']) and int(avals[0]) > 3):
         continue

      # Correct forward strand locus offset.
      if dirflag:
         for i,symb in enumerate(asymb):
            if symb == 'M' or symb == 'D':
               rloc += int(avals[i])

      # Add locus to list for current barcode
      if maps.has_key(barcode):
         if maps[barcode].has_key(rchr):
            maps[barcode][rchr].append(rloc if dirflag else -rloc)
         else:
            maps[barcode][rchr] = [rloc if dirflag else -rloc,]
      else:
         d = {}
         d[rchr] = [rloc if dirflag else -rloc,]
         maps[barcode] = d

#DEBUG
"""      
print "loci-barcode table:"
for loc in maps:
   sys.stdout.write(maps[loc][0])
   bcds = maps[loc][1]
   for bcd in bcds:
      sys.stdout.write('\t('+str(bcds[bcd])+')'+bcd)
   sys.stdout.write('\n')
"""

bcd_bestloc = {}
loc_bestbcd = {}

# Merge loci
for brcd in maps:
   bcd_best  = None
   best_cnt  = 0
   loc_cnt   = 0
   bcd_tot   = 0
   for c in maps[brcd]:
      loc_bcd = maps[brcd][c]
      loc_bcd = sorted(loc_bcd)
      # Add fake element to trigger push of last locus
      loc_bcd.append(loc_bcd[-1] + loc_merge_dist + 2)
      last = loc_bcd[0]
      best = loc_bcd[0]
      lcnt = 1
      bcnt = 1
      cnt  = 1
      for l in loc_bcd[1:]:
         if l - last == 0:
            lcnt += 1
         else:
            if lcnt > bcnt:
               best = last
               bcnt = lcnt
            if l - last <= loc_merge_dist:
               cnt += lcnt
               lcnt = 1
            else:
               # Final count of merged locus
               cnt += lcnt
               # 'best' contains canonical position
               loc_key = c+':'+str(abs(best))+(':+' if best > 0 else ':-')
               # Check whether this is the locus with more counts for current brcd
               if cnt > best_cnt:
                  bcd_best = loc_key
                  best_cnt = cnt
               # Update best barcode for current locus, format:
               # [0] best barcode
               # [1] reads of the best barcode
               # [2] total reads for this locus
               # [3] number of different barcodes
               if loc_bestbcd.has_key(loc_key):
                  loc_bestbcd[loc_key][2] += cnt
                  loc_bestbcd[loc_key][3] += 1
                  if cnt > loc_bestbcd[loc_key][1]:
                     loc_bestbcd[loc_key][0] = brcd
                     loc_bestbcd[loc_key][1] = cnt
               else:
                  loc_bestbcd[loc_key] = [brcd,cnt,cnt,1]

               # Update barcode reads
               loc_cnt += 1
               bcd_tot += cnt
               # Reset merge variables
               best = l
               bcnt = 1
               lcnt = 1
               cnt  = 1
         last = l
   # For each barcode store:
   # [0] best locus (more reads)
   # [1] reads of the best locus
   # [2] total reads of this barcode
   # [3] number of different loci
   bcd_bestloc[brcd] = [bcd_best, best_cnt, bcd_tot, loc_cnt]

del maps

bcdlist = sorted(bcd_bestloc.keys())

recomb_file.write('brcd\tchr\tloc\tstrand\treads\tbcd_tot\tbcd_numloc\tloc_tot\tloc_numbcd\tscore\n')

for brcd in bcdlist:
   # Recover barcode info
   bcd_match = bcd_bestloc[brcd]
   if bcd_match[0] is None:
      continue

   # Check whether best locus and best barcode match each other.
   loc_match = loc_bestbcd[bcd_match[0]]
   if loc_match[0] != brcd:
      continue
   bcd_cnt = bcd_match[1]
   bcd_tot = bcd_match[2]
   bcd_num = bcd_match[3]
   loc_cnt = loc_match[1]
   loc_tot = loc_match[2]
   loc_num = loc_match[3]
   score = max((1-bcd_cnt*1.0/bcd_tot),(1-loc_cnt*1.0/loc_tot))
   if score != 0:
      score = int(round(max(0,-10*math.log10(score))))
   else:
      score = 100
   if loc_cnt >= min_reads and score >= min_score:
      integ_file.write(brcd+'\t'+bcd_match[0].replace(':','\t')+'\t'+str(loc_cnt)+'\t'+str(score)+'\n')
      
   # Write assignment info in recomb_file
   recomb_file.write(brcd+'\t'+bcd_match[0].replace(':','\t')+'\t'+str(loc_cnt)+'\t'+str(bcd_tot)+'\t'+str(bcd_num)+'\t'+str(loc_tot)+'\t'+str(loc_num)+'\t'+str(score)+'\n')

recomb_file.write('\nall loci found:\n')
for loc in loc_bestbcd.keys():
   locus = loc.split(':')
   if locus[0] == 'HIV':
      continue
   recomb_file.write('NA\t'+loc.replace(':','\t')+'\t'+str(loc_bestbcd[loc][2])+'\tNA\n')
