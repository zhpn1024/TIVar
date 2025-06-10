__all__=[]
import torch
import numpy
import os

motiflen = 16
lhead, ltail = 9, 4
char_to_int = {'A':[1,0,0,0], 'T':[0,1,0,0], 'C':[0,0,1,0], 'G':[0,0,0,1]}
char_to_int['U'] = char_to_int['T']
ulim, dlim = 0.8, 0.2
dth, fth = 0.05, 2/3.0

def encode(s):
  seq_encoded=[]
  if len(s) != motiflen: return []
  for l in s:
    l = l.upper()
    if l not in char_to_int: return []
    seq_encoded += char_to_int[l]
  return seq_encoded

lim = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
mlen = len(lim)
nns = None
path = os.path.dirname(__file__)
path += '/cv2best'
#path = os.path.abspath(os.path.dirname(__file__)) #os.getcwd()

def TI_format(sq):
  if sq is None: return 'None'
  return '{}-{}-{}'.format(sq[0:lhead], sq[lhead:lhead+3], sq[lhead+3:])

def m6cmp(seqs):
  global nns
  if nns is None:
    nns = [torch.load('{}_r5.model'.format(path), weights_only=False), torch.load('{}_r4.model'.format(path), weights_only=False), torch.load('{}_r3.model'.format(path), weights_only=False),
           torch.load('{}_r2.model'.format(path), weights_only=False), torch.load('{}_r1.model'.format(path), weights_only=False), torch.load('{}_r05.model'.format(path), weights_only=False),
          ]
  dj, x = {}, []
  for j, s in enumerate(seqs):
    c = encode(s)
    if len(c) == 0: continue
    dj[j] = len(x)
    x.append(c)
  #x = [encode(s) for s in seqs] # [encode(s1), encode(s2)]
  xx = torch.tensor(x, dtype=torch.float, requires_grad=True)
  preds = [nns[i](xx) for i in range(mlen)]

  res = []
  for j in range(len(seqs)):
    est, outrange, used = {}, {}, {}
    ps = [None] * mlen
    if j not in dj:
      res.append([None, ps, used])
      continue
    j1 = dj[j]
    for i in range(mlen):
      p = float(preds[i][j1])
      ps[i] = p
      #print(lim[i], round(p, 4), round(p * lim[i], 4))
      est[i] = p * lim[i]
      if i != 0 and p > ulim: outrange[i] = 1
      elif i != 5 and p < dlim: outrange[i] = -1
      else: outrange[i] = 0
      if outrange[i] == 0: used[i] = est[i]

    #consist = False
    l = len(used)
    if l > 0:
      m = numpy.median(list(used.values()))
    else: m = None
    while l < mlen:
      for i in range(mlen):
        if i in used: continue
        if m is None: used[i] = est[i]
        elif outrange[i] == 1 and m <= ulim * lim[i]: used[i] = est[i]
        elif outrange[i] == -1 and m >= dlim * lim[i]: used[i] = est[i]
      m = numpy.median(list(used.values()))
      if len(used) == l: break
      l = len(used)
    #print(seq[j], round(m, 4), used)
    res.append([m, ps, used])
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
    else: s = 'No_TIS'
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

