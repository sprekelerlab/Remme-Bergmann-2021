function weights = bimodal_dist(mu, Ne_half)
    weights = [min(exprnd(mu, Ne_half, 1), 1); (1 - min(exprnd(mu, Ne_half, 1), 1) )];