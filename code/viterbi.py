import numpy as np

def viterbi(spec):
    """
    Run the Viterbi dynamic programming recursion for an HMM frequency-tracking problem.

    Parameters
    ----------
    spec : np.ndarray, shape (Nf, NT)
        "Emission" log-likelihoods (or scores) for each frequency bin f at each time step t.
        In this repo, spec is typically derived from a per-segment detection statistic
        (e.g., a periodogram / matched-filter score) mapped onto a frequency-bin grid.

        - Nf = number of frequency bins (states)
        - NT = number of time segments (time steps)

    Returns
    -------
    delta : np.ndarray, shape (Nf, NT)
        Dynamic programming table of best (max) cumulative score ending in state f at time t.
        delta[f, t] = best score of any path that ends at frequency bin f at time step t.

    backptrs : np.ndarray, shape (Nf, NT)
        Backpointer table storing the argmax predecessor state for each (f, t).
        backptrs[f, t] = previous frequency bin (state) at time t-1 that yields delta[f, t].
    """
    # Nf: number of frequency bins (states)
    # NT: number of time steps (segments)
    Nf, NT = spec.shape
    
    # delta holds the best cumulative score up to each (f, t)
    delta = np.zeros(spec.shape)

    # backptrs holds the predecessor state index to reconstruct best path later
    backptrs = np.zeros(spec.shape)

    # Initialisation at t=0:
    # Best path ending at frequency f at t=0 is just the local score at (f, 0)
    delta[:, 0] = spec[:, 0]

    # Recursion: fill delta and backptrs forward in time
    for t in range(1, NT):
        # Handle "interior" frequency bins that have 3 possible predecessors:
        # f-1 (down), f (stay), f+1 (up)
        # This encodes a transition model: the frequency can move by at most 1 bin per step.
        for f in range(1, Nf-1):
            # Consider the three candidate predecessor states at time t-1
            prev_window = delta[f - 1:f + 2, t - 1]  # [f-1, f, f+1] at previous time
            
            # Best cumulative score = best predecessor score + current emission score
            delta[f, t] = np.max(prev_window) + spec[f, t]

            # Store which predecessor state achieved the max:
            # argmax(prev_window) returns 0,1,2 => map back to global state index (f-1, f, f+1)
            backptrs[f, t] = f - 1 + np.argmax(prev_window)

        # Boundary case: f = 0 has only two predecessors (0 and 1),
        # because f-1 would be out of bounds.
        prev_window = delta[0:2, t-1]
        delta[0, t] = np.max(prev_window) + spec[0, t]
        backptrs[0, t] = np.argmax(prev_window)

        # Boundary case: f = Nf-1 has only two predecessors (Nf-2 and Nf-1),
        # because f+1 would be out of bounds.
        prev_window = delta[-2:, t-1]
        delta[-1, t] = np.max(prev_window) + spec[-1, t]
        backptrs[-1, t] = Nf - 2 + np.argmax(prev_window)

    # Note: this function does not pick the final best ending state.
    # Typically you choose start = argmax(delta[:, -1]) and then backtrace.
    return delta, backptrs

def backtrace(backptrs, start):
    """
    Reconstruct the most likely state (frequency-bin) path from backpointers.

    Parameters
    ----------
    backptrs : np.ndarray, shape (Nf, NT)
        Backpointer table from viterbi(), where backptrs[f, t] indicates the best predecessor
        state at time t-1 for a path ending at (f, t).

    start : int
        The final state (frequency bin) to start backtracing from at time t = NT-1.
        Common choice: start = np.argmax(delta[:, -1])

    Returns
    -------
    path : list[int], length NT
        The recovered best path of frequency-bin indices from t=0..NT-1.
        path[t] is the frequency bin at time step t.
    """
    # Initialise path with the final state at the final time step
    path = [int(start)]
    Nf, NT = backptrs.shape

    # Walk backwards from t = NT-1 down to 1.
    # At each step, prepend the predecessor state stored in backptrs.
    for t in range(NT - 1, 0, -1):
        prev_state = int(backptrs[path[0], t])  # predecessor of current first element at time t
        path.insert(0, prev_state)              # prepend to build path from t=0 upwards

    return path
