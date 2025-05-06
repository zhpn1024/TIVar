#!/bin/env python
from tivar.zbio import io, tools, mut, fa
from tivar import lib
import traceback

def set_parser(parser):
  parser.add_argument('-i', dest='input', required=True, type=str, help='input vcf file')
  parser.add_argument("-g", type=str, dest="genepath", required=True, help='Gene annotation file')
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output result file")
  parser.add_argument("-f", type=str, dest="genomefapath", required=True, help="Genome fasta file")
  parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if vcf chr ids are not same as chr ids in gtf/fasta files")

global chrmap
chrmap = {}

def vcfIter(fn, chrmap = chrmap):
  for lst in io.splitIter(fn):
    if lst[0].startswith('#'): continue
    if lst[4][0] not in ('A', 'T', 'C', 'G'):
      print('Skip unclear alt site {}'.format(lst[0:5]))
      continue
    if lst[4].find(',') > 0:
      print('Skip multiallele site {}'.format(lst[0:5]))
      continue
    m = mut.Mut(lst[0], int(lst[1])-1, lst[4], len(lst[3]))
    m.start, m.stop = m.pos, m.end
    m.tag = '{}:{}:{}>{}'.format(lst[0], lst[1], lst[3], lst[4])
    if m.chr in chrmap: m.chr = chrmap[m.chr]
    yield m


def run(args):

  #chrmap = {}
  if args.chrmap is not None :
    for lst in io.splitIter(args.chrmap, sep=None):
      if len(lst) < 2: continue
      chrmap[lst[0]] = lst[1]
      chrmap[lst[1]] = lst[0]
    fa.chrmap = chrmap

  genome = fa.Fa(args.genomefapath)
  mg = mut.MutGenome(args.genomefapath)

  #print('Loading gene annotation...')
  #anno = {}
  #for g in io.geneIter(args.genepath):
  #  if g.chr not in anno: anno[g.chr] = []
  #  anno[g.chr].append(g)
  #for chr in anno: pass

  #print('Processing...')

  motiflen = lib.motiflen # 16
  lhead, ltail = lib.lhead, lib.ltail # 9, 4

  outfile = open(args.output, 'w')
  outfile.write('Gid\tTid\tVar\tGenoPos\tStrand\tPos\tRefSeq\tAltSeq\tEffeRef\tEffeAlt\tDiff\tFC\tType\n')

  for m, g in tools.overlap_iter(vcfIter(args.input), io.geneIter(args.genepath)):
    #print(m.tag, g.id)
    mg.reset_mut()
    mg.add_mut(m)
    lindel = m.indel_len()
    used = {}
    #print(g.id, m)
    for t in g.trans:
      #print(t.id)
      tsq = genome.transSeq(t)
      msq = mg.transSeq(t)
      if tsq == msq: continue
      if t.strand != '-': p0 = t.cdna_pos(m.pos) # variant pos in trans
      else: p0 = t.cdna_pos(m.end) 
      if p0 is None: continue
      #if tsq[p0] == msq[p0]: p0 += 1
      sqref, sqalt, sqcmp = [], [], []
      p1, p2 = p0 - motiflen + 1, p0 + 1 # motif scan range
      if p1 < 0: p1 = 0
      if p2 + motiflen > len(tsq): p2 = len(tsq) - motiflen
      sqref = [tsq[i:i+motiflen] for i in range(p1, p2)]
      sqalt = [msq[i:i+motiflen] for i in range(p1, p2+lindel)]
      l1, l2 = len(sqref), len(sqalt)
      if l1 == 0 or l2 == 0: continue
      if lindel == 0:
        for i in range(0, l1):
          if len(sqref[i]) < motiflen or len(sqalt[i]) < motiflen: continue
          sqcmp.append([i, i])
      else:
        d = p0 - p1 - lhead
        if d > 0: # p0 - p1 > lhead:
          for i in range(0, min(d, l2)):
            if i >= l2 or len(sqalt[i]) < motiflen: continue
            sqcmp.append([i, i])
        if lindel < 0:
          if lindel < -2:
            d2 = abs(lindel) - 2
            for i in range(max(0, d), d+d2):
              if i >= l2 or len(sqalt[i]) < motiflen: continue
              if i == d and sqref[i][lhead] == sqalt[i][lhead]: sqcmp.append([i, i])
              else: sqcmp.append([i, None])
          for i in range(max(0, d+abs(lindel)-2), l1):
            if i >= l2 or len(sqalt[i]) < motiflen: continue
            sqcmp.append([i, i-abs(lindel)])
        elif lindel > 0:
          if lindel > 2:
            d2 = abs(lindel) - 2
            for i in range(max(0, d), d+d2):
              if i >= l1 or len(sqalt[i]) < motiflen: continue
              try:
                if i == d and sqref[i][lhead] == sqalt[i][lhead]: sqcmp.append([i, i])
                else: sqcmp.append([None, i])
              except Exception as e:
                print(e)
                print(i, l1, l2, sqref, sqalt)
                return False
          for i in range(max(0, d+abs(lindel)-2), l2): sqcmp.append([i-abs(lindel), i])

      seqs = sqref + sqalt
      for c in sqcmp:
        if c[1] is not None: c[1] += l1
      #print(seqs)
      try: res = lib.m6cmp(seqs)
      except Exception as e:
        print(e, str(e), type(e))
        traceback.print_exc()
        print(m.tag, g.id)
        print(l1, l2, seqs)
        return False
      #print(l1, l2, len(seqs), len(res), sqcmp)

      for c in sqcmp:
        r1, r2, sq1, sq2 = None, None, None, None
        is_na = False
        if c[0] is not None:
          r1, ps1, u1 = res[c[0]]
          sq1 = seqs[c[0]]
          pi = c[0]
          if res[c[0]][0] is None: is_na = True
        if c[1] is not None:
          r2, ps2, u2 = res[c[1]]
          sq2 = seqs[c[1]]
          if res[c[1]][0] is None: is_na = True
        if is_na: continue
        #r = [res[i] if i is not None else None for i in c]
        key = sq1, sq2
        if key not in used:
          ch, s, diff, fc = lib.getdiff(r1, r2)
          used[key] = ch, s, diff, fc
        else:
          ch, s, diff, fc = used[key]
          if not ch: continue
        if not ch: continue
        p = p1 + pi + lhead
        #print('{}\t{}\t{}\t{}\t{}'.format(t.gid, t.id, p, t.genome_pos(p), t.strand))
        #print(sq1)
        #print(sq2)
        #print(s)
        if r1 is not None: r1 = round(r1, 5)
        if r2 is not None: r2 = round(r2, 5)
        outfile.write(io.tabjoin(t.gid, t.id, m.tag, t.genome_pos(p), t.strand, p, lib.TI_format(sq1), lib.TI_format(sq2), r1, r2, round(diff,4), round(fc,4), s, '\n'))
  return True


if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())

