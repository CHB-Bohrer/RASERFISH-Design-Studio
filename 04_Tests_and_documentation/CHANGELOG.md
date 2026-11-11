Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

# 1.6.0 - 15 October 2026

- Verified ERRFI1 / DP1 startup, policy defaults and restore controls.
- Balanced 1 kb candidate sectors; 400 preferred + 200 extended candidates per region; bounded compatible-set repair.
- Standard compact results and standalone viewer; optional compressed full audit; portable hashed settings bundles.
- Scientific gates and exact comparison operators preserved from the saved ERRFI1 run.
- Reuse the installed v1.5 environment/databases while writing new runs in v1.6.
- 294 regression tests and 16 JavaScript pairing tests passed. Actual ERRFI1 generation produced 600 candidates, 30 per sector; full ERRFI1 QC/yield remains untested.

See RELEASE_NOTES.md and TESTING.txt for details and limits.

# 1.5.0 - 14 October 2026

Complete one-folder fresh installer; fixed-strand safeguard integrated; transcript-only RNA default with strict option; independent preset choices. See docs/history/v1_5/RELEASE_NOTES.md.

## 1.4.1 single-strand safeguard
- Require a fixed DNA-bound strand for the complete selection run; remove automatic per-gene choice.
- Retain RNA/complementarity exclusions and all non-strand thresholds.
- Recheck every selected site/repeat and mapped chromosome orientation before draft exports.
- Preserve old results with mixed-strand warnings; block new draft-order export for mixed/incomplete saved strand evidence.
- Add explicit strand-choice, preset-preservation, opposing-gene, repeat-copy, mapped-orientation, output-audit, and end-to-end tests.

# 1.4.0 - sequence-derived imager-binding graphics

- Added a live SVG showing actual antiparallel base pairing and 5-prime/3-prime ends.
- Gold is the imager; purple is the docking site; teal is the primary's target arm.
- Show two imager molecules on the two identical primary docking sites.
- Add a readout-preview selector and Show binding actions in Step 3.
- Reuse the diagram in the add-imager form and saved-candidate inspector.
- Mouse/keyboard position highlighting and self-contained SVG export.
- Keep both canonical stored/input sequences 5-prime to 3-prime.
- Retain all sequence, scientific-engine and policy inputs unchanged from 1.3.
- No new dependencies; no new experimental design or QC run.

# 1.3.0 - explicit imaging-oligo versus docking-site input

- Added a required input-type choice and both 5-prime-to-3-prime previews in step 3.
- Added an orientation review checkbox and an optional antiparallel display.
- Both front end and server resolve/validate the pair; legacy sequence fields keep
  their imager meaning. Conflicting fields are rejected, not overwritten.
- Detect duplicate pairs across input conventions; save input-type provenance.
- Show the pair in each saved candidate inspector; export explicitly named columns,
  imager_inventory.csv and readout_assignments.csv.
- Preserve existing target/oligo files and all scientific thresholds/defaults.

# 1.2.0 - 13 October 2026
- User-requested 45-55 C arm-Tm window with unchanged calculation assumptions.
- Adaptive short/extended length bands and independent, configurable QC budgets.
- Maximum compatible-site search compares baseline/infill/joint greedy attempts.
- Count-cap mode preserved; new length/selection evidence in UI, CSVs and JSON.
- Original paper preset, RNA safety, all other QC gates and dark interface retained.
- Original 163-probe order and target inputs unchanged.

# Design Studio changes

## 1.1.0 - 13 October 2026

- Initial arm range changed from 36-90 to paper-informed 36-42 nt. Default Tm window
  remains 42-50 C under the explicitly retained OligoMiner-style 2 nM/trace profile.
- Added field-level source provenance and an optional expanded G5-length preset.
  Primary GC limits are labeled project fallback; imager GC rules are not reused.
- Kept the user's five-base extended-footprint rule and equal adjacent bin planner.
- Added automatic gene-aware DNA-strand choice BEFORE candidate generation and
  checks for every repeat occurrence. Opposite/unknown gene directions exclude
  conflicting windows. Each selected primary exposes its actual bound strand.
- Synthetic references require user-confirmed transcription directions/intervals;
  missing direction is an error. Genomic input attempts gene/transcript retrieval.
- Added unspliced local gene RNA spans and annotated spliced RNA-feature screening.
  Exact complementary RNA tract threshold is 18 bases by default, independently of
  human DNA exemptions. RNA avoidance and transcript screening have separate scope.
- Added transcription/provenance outputs and explicit mixed-strand region warnings.
- Darkened the full Studio and saved results; preserved imager and annotation controls.
- Fixed candidate Tm-window display using each candidate's own Tm; fixed automatic
  orientation display so reverse-gene probes are labeled as targeting TOP (+).
- Replaced the preview's old demo with a clearly synthetic 3-bin RNA-direction
  tutorial. Did not change any of the original 163 project primary probes or targets.
- Added 45 tests. Total 120 tests: 118 passed, 2 real-engine tests skipped here.
  28 browser/API and 8 standalone checks passed, with documented sandbox browser transport bridge.

## 1.0.0 - earlier source version

QC-guided candidate designer, equal adjacent regions, imager dropdowns, intended
human-interval exemptions, feature adapters, and balanced greedy selection with
folding/pair and off-target gates. Retained unchanged in the previous delivered ZIP.
