__all__=[]

commands = {
  'predict': 'Predict potential TIS',
  'diff': 'Predict TI changes by VCF and transcript annotation',
  }

def load(cmd):
  if cmd == 'diff':
    from . import diff
    return diff
  elif cmd == 'predict':
    from . import predict
    return predict
  else:
    return None
