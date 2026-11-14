# geneweaver/cpu_baseline.py
"""Pure Python CRISPR off-target search (the CPU baseline).

A hit is a 20-base window that matches the guide with at most max_mm
mismatches and is followed by an NGG PAM (N = any base, G = code 2).
Genome and guide are uint8 arrays from parser.encode().
"""


def cpu_search(genome, guide, max_mm=3):
    """Return a list of (position, mismatches) for every hit in genome."""
    n, m = len(genome), len(guide)
    hits = []
    for i in range(n - m - 2):  # leave room for the 3-base PAM
        # PAM check first: it rejects ~15/16 of positions cheaply
        if genome[i + m + 1] != 2 or genome[i + m + 2] != 2:
            continue
        mm = 0
        for j in range(m):
            if genome[i + j] != guide[j]:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            hits.append((i, mm))
    return hits