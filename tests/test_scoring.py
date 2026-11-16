from geneweaver.scoring import off_target_score

GUIDE = "GACGCATAAAGATGAGACGC"


def mutate(guide, index):
    out = list(guide)
    out[index] = "A" if out[index] != "A" else "C"
    return "".join(out)


def test_exact_match_scores_100():
    assert off_target_score(GUIDE, GUIDE) == 100.0


def test_more_mismatches_lower_score():
    one = off_target_score(GUIDE, mutate(GUIDE, 10))
    two = off_target_score(GUIDE, mutate(mutate(GUIDE, 10), 15))
    assert two < one < 100


def test_mismatch_near_pam_scores_lower_than_far_from_pam():
    near_pam = off_target_score(GUIDE, mutate(GUIDE, 17))
    far_from_pam = off_target_score(GUIDE, mutate(GUIDE, 5))
    assert near_pam < far_from_pam