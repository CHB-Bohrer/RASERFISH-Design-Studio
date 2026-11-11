Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

SINGLE-STRAND SAFEGUARD — DESIGN STUDIO 1.6

The safeguard is already included. Do not apply an older strand patch.

Top (+) is the verified ERRFI1 default for ALL tracing regions. A fresh policy no
longer starts as "unselected". UI, backend, ERRFI1 preset and full ERRFI1 restore
agree. Bottom (-) remains an explicit choice for another design. Automatic,
mixed and unselected production modes remain invalid.

Applying ERRFI1 resets strand to Top (+) and the verified RNA policy while keeping
the current region/readout layout. Full ERRFI1 restore also returns the reference,
local 189019-209018 region and DP1. Imported settings retain their saved strand
until a reset is explicitly applied.

There is no per-gene switching. A top-targeting arm is the reverse complement of
the forward reference window. Exact repeat copies and every final selected site
must obey the trace-wide rule; mapped-fragment orientations are also audited.
Empty results are not a strand pass.

Transcript-sequence RNA screening is ON and blanket gene-span exclusion is OFF by
default. Strict span mode rejects conflicts instead of switching strands. Compact
output preserves the final strand audit.

Use Top (+) in every separate run for the current five-window ERRFI1 trace.
Separate runs still require combined-pool pair and strand QC. Sequence consistency
does not establish same-chromatid or maternal/paternal identity. This release has
not produced a newly screened full ERRFI1 probe set.
