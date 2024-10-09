__all__=[]
import torch
import numpy
import os

motiflen = 16
lhead, ltail = 9, 4
char_to_int = {'A':[1,0,0,0], 'T':[0,1,0,0], 'C':[0,0,1,0], 'G':[0,0,0,1]}
char_to_int['U'] = char_to_int['T']
ulim, dlim = 0.9, 0.1
dth, fth = 0.05, 2/3.0

def encode(s):
  seq_encoded=[]
  for l in s:
    seq_encoded += char_to_int[l.upper()]
  return seq_encoded

lim = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
nns = None
path = os.path.dirname(__file__)
#path = os.path.abspath(os.path.dirname(__file__)) #os.getcwd()

def TI_format(sq):
  if sq is None: return 'None'
  return '{}-{}-{}'.format(sq[0:lhead], sq[lhead:lhead+3], sq[lhead+3:])

def m6cmp(seqs):
  global nns
  if nns is None:
    nns = [torch.load('{}/my_nn_6_r5.model'.format(path)), torch.load('{}/my_nn_6_r4.model'.format(path)), torch.load('{}/my_nn_6_r3.model'.format(path)),
           torch.load('{}/my_nn_6_r2.model'.format(path)), torch.load('{}/my_nn_6_r1.model'.format(path)), torch.load('{}/my_nn_6_r05.model'.format(path)),
          ]
  x = [encode(s) for s in seqs] # [encode(s1), encode(s2)]
  xx = torch.tensor(x, dtype=torch.float, requires_grad=True)
  preds = [nns[i](xx) for i in range(6)]

  res = []
  ps = [None] * 6
  for j in range(len(seqs)):
    est, outrange, used = {}, {}, {}
    for i in range(6):
      p = float(preds[i][j])
      ps[i] = p
      #print(lim[i], round(p, 4), round(p * lim[i], 4))
      est[i] = p * lim[i]
      if i != 0 and p > ulim: outrange[i] = 1
      elif i != 5 and p < dlim: outrange[i] = -1
      else: outrange[i] = 0
      if outrange[i] == 0: used[i] = est[i]

    if len(used) > 0:
      m = numpy.median(list(used.values()))
    else: m = None
    #print(round(m, 4), used)
    for i in range(6):
      if i in used: continue
      if m is None: used[i] = est[i]
      elif outrange[i] == 1 and m <= ulim * lim[i]: used[i] = est[i]
      elif outrange[i] == -1 and m >= dlim * lim[i]: used[i] = est[i]
    m = numpy.median(list(used.values()))
    #print(seq[j], round(m, 4), used)
    res.append([m, ps])
  return res #, diff, fc

def getdiff(r1, r2, dth = dth, fth = fth):
  ch, diff, fc = False, None, 0
  if r1 is None: # and r2 > dth:
    #s = 'T1=None, T2={0}, Diff={0}, FC=0'.format(round(r2, 4))
    diff = r2
    if r2 > dth: s = 'TIS_emerge'; ch = True
    else: s = 'No_TIS'
  elif r2 is None: # and r1 > dth:
    #s = 'T1={0}, T2=None, Diff={0}, FC=0, TIS loss'.format(round(r1, 4))
    diff = -r1
    if r1 > dth: s = 'TIS_loss'; ch = True
    else: s += 'No_TIS'
  else:
    diff = r2 - r1
    if diff == 0: fc = 1
    elif diff > 0: fc = r1 / r2
    else: fc = r2 / r1
    #s = 'T1={}, T2={}, Diff={}, FC={}'.format(round(r1, 4), round(r2, 4), round(diff, 4), round(fc, 4))
    if abs(diff) < dth or fc > fth: s = 'TI_no_change'
    elif diff > 0: s = 'TI_increased'; ch = True
    else: s = 'TI_decreased'.format(diff, fc); ch = True
  return ch, s, diff, fc

