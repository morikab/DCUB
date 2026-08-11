import numpy as np

from modules.sequence_family.clustering_optimization import create_n_clusters, make_distance_matrix

# Fixed, hand-picked 4-organism x 5-codon CAI-profile-shaped matrix: two
# clearly-similar pairs of rows, so clustering into 2 groups has an
# unambiguous expected answer regardless of scikit-learn's internal tie-
# breaking. Values and expected outputs below were captured by running this
# exact input against the current (pre-migration) code.
CLUSTERING_MAT = np.array([
    [0.9, 0.1, 0.8, 0.2, 0.7],
    [0.85, 0.15, 0.75, 0.25, 0.65],
    [0.1, 0.9, 0.2, 0.8, 0.3],
    [0.15, 0.85, 0.25, 0.75, 0.35],
])


def test_make_distance_matrix_is_unchanged():
    distance_matrix = make_distance_matrix(CLUSTERING_MAT)
    # Expected near-zero cells (1.11022302e-16) are float64 rounding noise from
    # spearmanr's internal Pearson-on-ranks computation; the true mathematical
    # value is 0.0 (perfectly correlated rows), but summation order in floating-
    # point arithmetic produces this tiny epsilon. Use atol=1e-9 to tolerate
    # rounding-noise drift across numpy/scipy versions while still catching real
    # behavior changes (actual dissimilarity is 2.0, far larger than the tolerance).
    expected = np.array([
        [1.11022302e-16, 1.11022302e-16, 2.0, 2.0],
        [1.11022302e-16, 1.11022302e-16, 2.0, 2.0],
        [2.0, 2.0, 1.11022302e-16, 1.11022302e-16],
        [2.0, 2.0, 1.11022302e-16, 1.11022302e-16],
    ])
    np.testing.assert_allclose(distance_matrix, expected, atol=1e-9)


def test_create_n_clusters_groups_similar_organisms_together():
    labels = create_n_clusters(CLUSTERING_MAT, n_clus=2)
    # Rows 0,1 (the "0.9/0.1..." pair) must land in one cluster; rows 2,3
    # (the "0.1/0.9..." pair) must land in the other. Cluster label
    # *numbering* is arbitrary (scikit-learn doesn't guarantee 0 vs 1 means
    # anything) — only the grouping matters.
    assert labels[0] == labels[1]
    assert labels[2] == labels[3]
    assert labels[0] != labels[2]
