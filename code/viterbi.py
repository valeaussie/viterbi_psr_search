import numpy as np

def viterbi(spec):
    Nf, NT = spec.shape
    
    delta = np.zeros(spec.shape)
    backptrs = np.zeros(spec.shape)

    delta[:, 0] = spec[:, 0]

    for t in range(1, NT):
        for f in range(1, Nf-1):
            delta[f, t] = np.max(delta[f-1:f+2, t-1]) + spec[f, t]
            backptrs[f, t] = f - 1 + np.argmax(delta[f-1:f+2, t-1])

        delta[0, t] = np.max(delta[0:2, t-1]) + spec[0, t]
        backptrs[0, t] = np.argmax(delta[0:2, t-1])

        delta[-1, t] = np.max(delta[-2:, t-1]) + spec[-1, t]
        backptrs[-1, t] = Nf - 2 + np.argmax(delta[-2:, t-1])
    
    return delta, backptrs

def backtrace(backptrs, start):
    path = [int(start)]
    Nf, NT = backptrs.shape
    for t in range(NT-1, 0, -1):
        path.insert(0,int(backptrs[path[0], t]))

    return path
