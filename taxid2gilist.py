#! /usr/bin/env python
from __future__ import division, print_function
import sys
from Bio import Entrez
Entrez.email = 'harald.detering@fu-berlin.de'

# make sure required parameters have been provided
if len(sys.argv) < 2:
  print("usage: %s taxid" % sys.argv[0], file=sys.stderr)
  sys.exit(1)

# get taxonomy ID from command line params
taxid = sys.argv[1]

# query NCBI
query = "txid%s[Organism:exp]" % taxid
handle = Entrez.esearch(db='nucleotide', term=query, usehistory=True)
search_result = Entrez.read(handle)
webenv = search_result['WebEnv']
query_key = search_result['QueryKey']
record_cnt = int(search_result['Count'])
print('Found %s records for taxid %s.' % (record_cnt, taxid), file=sys.stderr)

# get records
batch_size = 20000
print('Downloading records:', file=sys.stderr)
for start in range(0, record_cnt, batch_size):
  print('\trecord %s to %s' % (start, start+batch_size-1), file=sys.stderr)
  handle = Entrez.esummary(db='nucleotide', webenv=webenv, query_key=query_key, 
                           retstart=start, retmax=batch_size)
  for record in Entrez.parse(handle):
    print('%s' % record['Gi'])
  
