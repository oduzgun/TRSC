# Post-treatment implementation boundary

The accepted identification and inference layers are implemented in package
R/identification.R, R/certificates.R and R/inference.R. See the public help and
executable vignettes. This directory is reserved for future research artifacts,
not an alternative estimator.

Population sets, plug-in sets and confidence regions have distinct labels.
Causal bridges, continuation orders/radii, auxiliary exclusions and statistical
certificates remain explicit supplied conditions. No automatic treatment-effect
identification or post-treatment treated-outcome fitting is performed.
