"""Off-target scoring (MIT specificity score, Hsu et al. 2013).

Higher score = the target is more likely to be cut (a more dangerous
off-target). Mismatches near the PAM (the seed region) lower the score a
lot because Cas9 tolerates them poorly; mismatches far from the PAM are
tolerated, so they lower the score less.
"""

# Mismatch weights for guide positions 1..20 (position 20 touches the PAM).
W = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
     0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]


def off_target_score(guide: str, target: str) -> float:
    """Return a 0-100 score for a 20-base target against a 20-base guide."""
    mism = [i for i in range(20) if guide[i] != target[i]]
    if not mism:
        return 100.0
    s1 = 1.0
    for i in mism:
        s1 *= 1 - W[i]
    n = len(mism)
    d = 19 if n == 1 else (mism[-1] - mism[0]) / (n - 1)
    s2 = 1 / (((19 - d) / 19) * 4 + 1)
    s3 = 1 / (n * n) if n > 1 else 1
    return 100 * s1 * s2 * s3
