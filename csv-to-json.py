#!/usr/bin/env python3

"""
Convert leaf resolution data from csv to JSON.
input file format: name, resolution

output: JSON dictionary mapping name to a dictionary of parameters
  where the only parameter is px_mm

If the output file exists, it is read in initially, so that existing data is preserved.
Existing leaves are updated, with a message printed if the resolution changed.

Usage: 
python3 csv-to-json.py [input.csv] [output.json]

"""

import csv
import json
import sys

# take input file name from the command line if given
if len(sys.argv) > 1:
    inFileName = sys.argv[1]
else:
    inFileName = 'resolutions.csv'

# take output file name from the command line if given
if len(sys.argv) > 2:
    outFileName = sys.argv[2]
else:
    outFileName = 'leaves.json'

# read the database if it exists
fileHandle = None
leaves = {}
try:
    fileHandle = open(outFileName, 'r')
    print('Reading existing leaf data from ' + outFileName + '.')
except:
    print('Could not load the leaf data base from ' + outFileName + '. Creating an empty data base.')

if fileHandle != None:
    # if the file exists, but the JSON is not correct, this will fail
    # This is on purpose, since if we continue and save the database at the end,
    # the data in it will be overwritten
    leaves = json.load(fileHandle)
    

incsv = csv.reader(open(inFileName, "r"))
for row in incsv:
    if row[0] != '' and row[1] != '':
        name = row[0]
        px_mm = float(row[1])

        if name not in leaves:
            leaves[name] = {}
            
        if 'px_mm' in leaves[name]:
            if leaves[name]['px_mm'] != px_mm :
                print ('Updating resolution of ' + name + "  old:%.1f  new:%.1f"%(leaves[name]['px_mm'], px_mm))

        leaves[name]['px_mm'] = px_mm

#    else:
##      An incomplete row in the input csv. We ignore it. 
#       print ('Skipping row:'+ str(row))

of = open(outFileName, "wt")
json.dump (leaves, of, sort_keys=True, indent=2, separators=(',', ': '))
