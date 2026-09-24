# Package completion plan

1. Audit and record the frozen-paper crosswalk; reproduce baseline tests.
2. Test and repair metric-factor consistency. Add validated conic problem
   assembly, exact-arithmetic upper certification and status diagnostics.
3. Implement evaluation fibers, anchored dynamic sets/support/minimum radius,
   sensitivity and same-witness auxiliary refinement; verify analytic cases.
4. Implement explicit anchor/geometry certificates, joint q2 uncertainty,
   outer support, selection probability accounting, effect functionals and
   deployment modulus; add theorem-aligned and adversarial tests.
5. Document all new interfaces; retain the three existing vignettes and add
   four focused executable guides and two post-treatment examples.
6. Extend the research harness with separated observable/truth inputs and
   post-treatment discrimination cases; freeze SIMULATION_DESIGN before any
   large Monte Carlo execution. Run small smoke tests only.
7. Regenerate help, update CI, execute all tests/examples/vignettes, build and
   check the package. Inspect exports and final diff; document limitations.
8. Publish disciplined commits on `astra/package-completion` and open one PR
   into main. Do not merge.

API design: preserve existing pre-treatment calls. Represent a post-treatment
set as a lazy convex support oracle with explicit geometry, named deployment
matrix, anchor representatives and structural q/R. Anchor representatives
fix profiles only. Inference replaces fixed anchors with joint feasible lifts.
Auxiliary maps act on those same future witnesses. Never silently select a
common lift or remove prediction-null directions. Certificate constructors
require finite values and provenance; they do not manufacture population
bounds from empirical diagnostics.

Implementation is authorized by the user's complete specification. Routine
architecture decisions and reversible edits proceed without another approval
round. Meaningful tests are written and observed failing before new code.
