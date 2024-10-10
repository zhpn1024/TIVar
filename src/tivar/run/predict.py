#!/bin/env python
from tivar.zbio import io, tools, mut, fa
from tivar import lib

def set_parser(parser):
  parser.add_argument('-s', dest='seq', type=str, help='input squence')
  parser.add_argument('-S', dest='seqfile', type=str, help='input squence text file')
  parser.add_argument("-g", type=str, dest="genepath", help='Gene annotation file')
  parser.add_argument("-f", type=str, dest="genomefapath", help="Genome fasta file")
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output result file")
  parser.add_argument("--th", type=str, help="Output score threshold, default 0.05")
  parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if vcf chr ids are not same as chr ids in gtf/fasta files")

global chrmap
chrmap = {}

def seqIter(seq):
  yield 'Seq', seq

def seqFileIter(seqfile):
  id, sq = 'Seq', ''
  for l in open(seqfile):
    l = l.strip()
    if l == '' : continue
    if l.startswith('>'):
      if sq != '': yield (id, sq)
      sq = ''
      id = l[1:]
    else:
      l = l.replace('U','T').replace('u','t')
      sq += l
  yield (id, sq)

def transIter(genomefapath, genepath):
  genome = fa.Fa(genomefapath)
  for g in io.geneIter(genepath):
    for t in g.trans:
      sq = genome.transSeq(t)
      id = '{} {}'.format(t.id, g.id)
      yield (id, sq)

def vcfIter(fn, chrmap = chrmap):
  for lst in io.splitIter(fn):
    if lst[0].startswith('#'): continue
    if lst[4].find(',') > 0:
      print('Skip multiallele site {}'.format(lst[0:5]))
      continue
    m = mut.Mut(lst[0], int(lst[1])-1, lst[4], len(lst[3]))
    m.start, m.stop = m.pos, m.end
    if m.chr in chrmap: m.chr = chrmap[m.chr]
    yield m


def run(args):

  if args.seq is not None: sqiter = seqIter(args.seq)
  elif args.seqfile is not None: sqiter = seqFileIter(args.seqfile)
  elif args.genomefapath is not None and args.genepath is not None: sqiter = transIter(args.genomefapath, args.genepath)
  else:
    print('Invalid input! Please provide either -s (squence only), -S (sequence file) or (-g annotation -f genome) combined input.')
    exit(-1)

  global chrmap
  if args.chrmap is not None :
    for lst in io.splitIter(args.chrmap, sep=None):
      if len(lst) < 2: continue
      chrmap[lst[0]] = lst[1]
      chrmap[lst[1]] = lst[0]
    fa.chrmap = chrmap

  motiflen = lib.motiflen
  lhead, ltail = lib.lhead, lib.ltail
  dth = lib.dth
  if args.th is not None: dth = args.th

  outfile = open(args.output, 'w')
  outfile.write('SeqID\tPos\tStartSeq\tEffScore\n')
  for id, sq in sqiter:
    print('Predicting {}...'.format(id))
    p1, p2 = 0, len(sq) - motiflen
    if p2 < 1: continue
    msqs = [sq[i:i+motiflen] for i in range(p1, p2)]
    res = lib.m6cmp(msqs)
    for i, r in enumerate(res):
      m = r[0]
      if m >= dth:
        msq = msqs[i]
        sqstr = '{}-{}-{}'.format(msq[0:lhead], msq[lhead:lhead+3], msq[lhead+3:])
        outfile.write(io.tabjoin(id, i+lhead, sqstr, round(m, 5), '\n'))
  return True


if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())

