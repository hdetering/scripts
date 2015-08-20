#!/usr/bin/env python3
'''
author        | Harald Detering (harald.detering@gmail.com)
last modified | 2015-08-19
'''

import argparse, gzip, pysam, sys

def print_filtered_reads(reads, alignment, refs, exclusive):
  read_count = 0
  for m in alignment:
    if (m.rname in refs or m.mrnm in refs) and exclusive:
      continue
    elif not exclusive:
      continue
    
    if not m.is_supplementary:
      read_count += 1
      print("@%s/%s\n%s\n+\n%s" % (m.qname, 2 if m.is_read2 else 1, m.query, m.qqual))
  
  return read_count

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-r', '--reads', required=True, help="FASTQ file containing reads (may be gzip'ed)")
  parser.add_argument('-m', '--mapping', required=True, help="SAM/BAM file with read mappings")
  parser.add_argument('-e', '--exclude', action='store_true', help="exclusive filter: reads mapping to specified references are excluded in the output (default: inclusive)")
  parser.add_argument('-k', '--keys', help="reads mapping to these references (comma-separated) are included/excluded in the output")
  
  args = parser.parse_args()
  
  return args

def main(args):
  # process input params
  reads = gzip.open(args.reads)
  aln   = pysam.AlignmentFile(args.mapping, 'rb')
  # convert reference names to indices used by pysam
  ref_names = [x.strip() for x in args.keys.split(',')]
  ref_idx   = list(filter(lambda x: aln.references[x] in ref_names, range(len(ref_names))))
  
  print("Filtering reads in \n\t%s\n%smapping to references\n\t%s\nin alignment\n\t%s\n" % (args.reads, 'not ' if args.exclude else '', ','.join(ref_names), args.mapping), file=sys.stderr)
  read_cnt = print_filtered_reads(reads, aln, ref_idx, args.exclude)
  print("Written %s reads to output" % read_cnt, file=sys.stderr)
  
  return True

if __name__ == '__main__':
  args = parse_args()
  main(args)
