# Experimental post-treatment research boundary

This directory is reserved for post-treatment identification research that is not part of the stable TRSC API.

The current theory separates three questions:

1. **Predictive continuation:** what restrictions determine or bound a future predictive profile path?
2. **Evaluation ambiguity:** for a fixed profile, which values are compatible with different simplex lifts at an observed donor state?
3. **Cross-date lift structure:** what, if anything, links equivalent donor lifts across post-treatment dates?

For a predictive profile represented by `beta = B w`, the compatible lift set is

```text
W(beta) = { w in Delta_J : B w = beta }.
```

At a post-treatment donor vector `x`, the compatible evaluation set is

```text
Y(x; beta) = { x' w : w in W(beta) }.
```

This set is an interval whose endpoints are linear programs. These statements motivate future experimental evaluation-set code, but they do not by themselves identify a causal counterfactual.

Do not introduce stable `trsc_counterfactual()`, `trsc_treatment_effect()` or `trsc_att()` functions here until the continuation and causal-bridge assumptions are frozen.

Likewise, identifying or bounding an untreated conditional mean does not identify the realized untreated innovation. Empirical diagnostics must not be described as theorem verification.
