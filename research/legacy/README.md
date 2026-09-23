# Legacy TVSC implementation archive

These files are preserved to keep the software history auditable while the project migrates from TVSC to TRSC.

They are **not loaded by the TRSC package**.

## Classification

- `tvsc.blocks.R`, `tvsc.scale.R`, `tvsc.metrics.R`: version-1 infrastructure tied to the earlier preparation contract. Preserved as historical implementation evidence.
- `tvsc.vmetrics.R`: version-2 prefix scaling and predictor-metric construction. Useful implementation ideas remain, but its central predictor-metric abstraction is not the final representation-invariant predictive geometry.
- matching standalone tests: preserved with their original code and names.
- workflow-paper and workflow-strategy documents: preserved as historical design records, not current package documentation.

No legacy file should be sourced by stable TRSC code or CI except when a historical regression is being investigated explicitly.
