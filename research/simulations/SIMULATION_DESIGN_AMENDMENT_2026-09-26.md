# Simulation design operational amendment — 2026-09-26

## Status

This amendment was frozen before any large-run Monte Carlo performance findings were generated or inspected. It operationalizes choices that `SIMULATION_DESIGN.md` version 1 declared only partially. It does **not** revise TRSC theory, the package estimands, the California application, or the previously declared DGP axes.

The authoritative GitHub byte representation of `SIMULATION_DESIGN.md` has SHA-256 `7ba2f3086a1649751732a6853c7fff06825608ea6cc0435a6ab6397469b1cf51`. The historical `5a323b781a915ee313d5187c7d59d53a650947b4725f9af89b81ce42f678f6b8` refers to scientifically identical text without the terminal LF byte.

## Newly frozen operational specifications

The nested predictive libraries are frozen in `FROZEN_CANDIDATE_LIBRARIES.csv`: L1, L4, L12, and L24 exactly as specified in the Stage-A handoff. They are nested-library-complexity comparisons, not pure candidate-count experiments. Structural q remains separate. `FROZEN_DEPLOYMENT_EPSILON.csv` freezes epsilon/D at {0,0.01,0.1,0.5,1}; it is a reporting grid, not a DGP axis.

## Frozen registries

- `FROZEN_CELL_REGISTRY.csv`: SHA-256 `b69df091e331c6d93e786bf0e6300babb4a5d4b76e5cf6f09c3c242d16a51480`; 179 cells; 73,400 attempts.
- `FROZEN_CANDIDATE_LIBRARIES.csv`: SHA-256 `8f0a00b3ce3aef2ac294ab9946f2701f7977730ab9fa09038df0cb51e80f1be5`; 41 ordered membership rows.
- `FROZEN_DEPLOYMENT_EPSILON.csv`: SHA-256 `e5513f30939373cb98edfff718d2e74805095d97635741f407baf3abd029c017`; 5 reporting multipliers.

| Panel | Cells | Attempts |
|---|---:|---:|
| E0 | 13 | 9,000 |
| E1 | 12 | 12,000 |
| P0 | 36 | 11,200 |
| P1 | 10 | 5,200 |
| P2 | 3 | 3,000 |
| P3 | 35 | 7,000 |
| P4 | 44 | 8,800 |
| P5 | 11 | 2,200 |
| P6 | 8 | 8,000 |
| P7 | 3 | 3,000 |
| P8 | 4 | 4,000 |
| **Total** | **179** | **73,400** |

The seed rule remains `1000000 + 10000 * scenario_index + repetition`. Equal seeds across sensitivity cells are common-seed coupling only; CRN language requires demonstrably identical primitive draws. Attempt identity is `registry_hash + cell_id + repetition`.

Requested and effective dimensions are recorded separately; runtime execution must verify generator metadata. Population Rmin remains explicitly design/oracle information and is never represented as feasible data-driven tuning.

## Separate gates

The executable registry gate and the nontrivial nominal-inference gate are separate. A reproducible estimator/identification/oracle/plugin/conservative-diameter study may be execution-ready while the theorem-based nominal risk protocol remains NOT READY because numerical theorem inputs are unavailable.

## Immutability

After scientific execution begins, these registries are immutable. Any later change requires a dated amendment retaining the original registries and all prior outputs. Unsupported and failed cells remain in attempted-repetition accounting.


## Pre-results effective-rank correction

Before scientific execution, the registry's `rank_effective` field was audited against the generated population predictive geometry rather than the latent-dimension input. The post-treatment generators `exact_ambiguity`, `near_redundancy`, and `weak_excitation` with requested rank 3 intentionally overwrite one factor-loading direction, so their actual population predictive rank is 2. The 15 affected P0/P2 rows now record `rank_requested=3`, `rank_effective=2`, and an explicit generator-override note. This changes registry metadata and its hash only; it does not change generator code, DGP values, seeds, cells, repetitions, or any result. No scientific repetitions had been executed when this correction was frozen.
