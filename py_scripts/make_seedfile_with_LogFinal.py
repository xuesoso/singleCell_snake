#!/oak/stanford/groups/quake/yuanxue/resources/anaconda3/bin/python

import sys
import fnmatch
import os

dir = sys.argv[1] # directory containing aligned samples (Log.final.out)
out = sys.argv[2]
if os.path.isdir(out) == False:
    out = os.path.dirname(os.path.abspath(out)) + '/seedfile.txt'
else:
    out = out + '/seedfile.txt'

matches = []
for root, dirnames, filenames in os.walk(dir):
    for filename in filenames:
        if 'Log.final.out' in filename:
            matches.append(os.path.abspath(root)) # directory of match

# Print directory containing each file
if os.path.isfile(out):
    os.remove(out)

with open(out, 'w') as f:
    for x in sorted(matches):
        f.write(x+'\n')
