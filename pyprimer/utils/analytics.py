import numpy as np
import tqdm
import matplotlib.pyplot as plt

def entropy_calculation(sequences, window_size=3, start=0, end=None, verbose=False, draw=False):
    if end is None:
        end = len(sequences['Sense Sequence'].values[0])
    
    length = end-start
    
    res = np.zeros(length)
    entropy = np.zeros(length)

    if not verbose:
        iterable = range(start, end)
    else:
        iterable = tqdm.trange(start, end)

    for i, nucl in enumerate(iterable):
        vals = []
        for sequence in sequences['Sense Sequence'].values:
            fragment = sequence[nucl:nucl+window_size]
            if 'N' not in fragment:
                vals.append(fragment)
        
        if len(set(vals))>1:
            res[i] = len(set(vals))
        
        ent = 0
        for val in list(set(vals)):
            p = vals.count(val)/len(vals)
            ent -= p*np.log(p)
        entropy[i] = ent

    return res, entropy