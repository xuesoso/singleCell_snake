#!/oak/stanford/groups/quake/yuanxue/resources/anaconda3/bin/python

import sys
import os

infile = sys.argv[1] # directory containing aligned samples (Fastqs)
out = sys.argv[2]
if os.path.isdir(out) == False:
    out = os.path.dirname(os.path.abspath(out)) + '/seedfile.txt'
else:
    out = out + '/seedfile.txt'
matches = []
for root, dirnames, filenames in os.walk(infile):
    r1, r2 = False, False
    for filename in filenames:
        if 'R1_001.fastq.gz' in filename or 'R1.fastq.gz' in filename:
            r1 = True
        elif 'R2_001.fastq.gz' in filename or 'R2.fastq.gz' in filename:
            r2 = True
    if r1 == True and r2 == True:
        matches.append(os.path.abspath(root)) # directory of match

# Print directory containing each file
if os.path.isfile(out):
    os.remove(out)

with open(out, 'w') as f:
    for x in sorted(matches):
        f.write(x+'\n')
