Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

> Historical feature notes. For the CURRENT 1.5 defaults, installation and RNA/strand policy, see START_HERE.txt and METHODS_AND_LIMITATIONS.md. Strict RNA exclusion is no longer the default.

# Design Studio 1.2: adaptive length and 45-55 C

## Requested change

The default on-target DNA-binding-arm Tm window is now strictly 45 < Tm < 55 C.
This is a USER setting, not a claim that the Ellenberg paper used this window.
The default Tm model remains OligoMiner-style DNA_NN3 at 390 mM sodium, 2 nM
probe / trace target, salt correction 5 and the 50% formamide x 0.65 C correction.
No incubation or folding temperature is changed by changing the selection window.

## Candidate generation

Preferred arms are 36-42 nt. Adaptive mode additionally examines 43-90 nt arms,
one nucleotide at a time, bounded by the region edge and adjustable arm ceiling.
The 90-nt ceiling is retained from the earlier broader project search; it is NOT
an optimal length or a directive to lengthen every oligo. Handles and spacers are
not included in these arm lengths.

Each candidate length must pass the configured GC, identical-base and dinucleotide
filters, the current Tm screen and annotated/declared RNA-direction checks. After
those, exact-copy and footprint tests are applied. At each start, each length band
retains up to two shortest locally compatible alternatives by default. Alternatives
failing the local tests do not consume those slots. Longer variants are useful when
short ones cannot reach the Tm window, or when a different boundary changes a QC
interaction. We do not assume increasing length resolves every off-target.

Preferred and extended candidates have separate spatially distributed shortlist
budgets: 160 and 80 per region initially. Extended candidates do not replace the
preferred shortlist. All retained candidates receive the SAME remaining local,
genomic/RNA, monomer, self-pair and selection-time pair checks. Full mode requires
the real engines/databases. Preview mode remains explicitly local-only.

These sampling limits mean not every possible length receives expensive QC. A
shortlisted longer arm may fail even though some unshortlisted alternative would
pass. Increase alternative/shortlist budgets to search more, rather than assuming
an underfilled bin has no possible additional probes.

## Selecting more compatible sites

Default maximum-site mode compares four actual feasible solutions:
1. Preferred-length-only baseline using packing order.
2. That baseline retained, with compatible additions from the entire candidate set.
3. Joint site-packing search over both length bands.
4. Joint spatial-balance search over both length bands.

At every addition, all exact repeat footprints are checked across all regions and
every already selected primary is checked for pair compatibility. The winner has
the largest exact intended-site count among these attempts. Ties prefer more filled
regions, a higher lowest-region count, more distinct primaries, spatial spread and
less total arm footprint, in that order. All selected pairs/sites are rechecked.
A preferred-only feasible baseline is included in the comparison, so the adaptive
winner cannot have fewer total sites than THAT SAME screened baseline.

This is a deterministic multi-start GREEDY heuristic. It is not a global optimum,
a statement that every design grows, or a prediction of experimental signal. A
new Tm window is not necessarily a superset of the old window: candidates below
45 C now fail even if they previously met 42-50 C.

The per-region target number is a benchmark, NOT a stop condition, in maximum-site
mode. The separately selectable capped target-count mode still stops at that count.
Region widths and boundaries do not grow to fit long candidates.

## What is unchanged

RNA avoidance and required direction confirmations; human intended-interval policy;
local/cross-bin checks; high/moderate DNA off-target cutoffs (42 / 35 C); RNA flags;
full-primary folding and pair rules; 42 C aqueous folding; five-base EXTENDED
footprint spacing; docking-TT-arm-TT-docking architecture; user imager assignments;
fixed equal adjacent tracing regions; dark interface and feature/annotation tools.
Longer candidates are independently re-evaluated, not merely relabeled old arms.
The existing 163 order sequences and original target references are not modified.

## UI and outputs

Step 4 exposes length strategy, preferred maximum, extension ceiling, additional
QC budget and count goal. Step 6 reports extended-primary counts, arm-length range,
per-candidate extension metadata and all compared selection attempts. New CSV
columns and selection_search.csv preserve these decisions. Original paper-informed
and earlier expanded-range presets remain available and labeled separately.

## Demonstration scope

The bundled synthetic tutorial uses 3 adjacent 1,000-base regions with constructed
forward/reverse/overlapping gene features. It is not a real human locus. Local-only
QC selected 33 primaries at 33 sites versus 25 sites for the preferred-only baseline;
15 selected primaries are longer than 42 nt. These counts demonstrate the algorithm,
not a predicted gain in the user's CFTR pool. No new full-reference CFTR design ran.
