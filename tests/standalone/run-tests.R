# Run from the repository root; these contracts intentionally use base R only.
source("tests/standalone/test-trsc-data.R", local = new.env())
source("tests/standalone/test-trsc-split.R", local = new.env())
cat("Standalone TRSC contracts passed.\n")
