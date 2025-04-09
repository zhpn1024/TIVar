#!/bin/env python
from tivar.zbio import io, tools, mut, fa
from tivar import lib

def set_parser(parser):
  parser.add_argument('-s', dest='seq', type=str, help='input squence')
  parser.add_argument('-S', dest='seqfile', type=str, help='input squence text file')
  parser.add_argument("-g", type=str, dest="genepath", help='Gene annotation file')
  parser.add_argument("-f", type=str, dest="genomefapath", help="Genome fasta file")
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output result file")
  parser.add_argument("-m", type=str, dest="model", help="Input model prefix")
  parser.add_argument("--th", type=float, help="Output score threshold, default 0.05")
  parser.add_argument("--all", action="store_true", help="Output all scores from 6 models")
  parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase output verbosity")
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
  if args.model is not None: lib.path = args.model

  outfile = open(args.output, 'w')
  if args.all: outfile.write('SeqID\tPos\tStartSeq\tEffScore\tM5\tM4\tM3\tM2\tM1\tM05\n')
  else: outfile.write('SeqID\tPos\tStartSeq\tEffScore\n')
  for id, sq in sqiter:
    if args.verbose > 0: print('Predicting {}...'.format(id))
    p1, p2 = 0, len(sq) - motiflen + 1
    if p2 < 1: continue
    msqs = [sq[i:i+motiflen] for i in range(p1, p2)]
    res = lib.m6cmp(msqs)
    for i, r in enumerate(res):
      m = r[0]
      if m >= dth:
        msq = msqs[i]
        sqstr = '{}-{}-{}'.format(msq[0:lhead], msq[lhead:lhead+3], msq[lhead+3:])
        if args.all:
          rstr = [round(r[1][j]*lib.lim[j], 5) for j in range(lib.mlen)]
          outfile.write(io.tabjoin(id, i+lhead, sqstr, round(m, 5), rstr, '\n'))
        else:
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

