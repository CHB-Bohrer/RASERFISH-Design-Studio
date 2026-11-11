Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

> Historical feature notes. For the CURRENT 1.5 defaults, installation and RNA/strand policy, see START_HERE.txt and METHODS_AND_LIMITATIONS.md. Strict RNA exclusion is no longer the default.

> HISTORICAL 1.1 CHANGE RECORD. The active 1.2 default supersedes the fixed
> length/Tm defaults below. See ADAPTIVE_LENGTH_CHANGES.md. RNA-avoidance rules
> remain unchanged. This file preserves the source explanation for the old preset.

# What changed in Design Studio 1.1

## What the initial preset means

The sources are separated at field level in `design_presets.py` and in step 4 of
the interface. Choosing a named preset is not a claim of reproducing every part of
a published pipeline. A saved result reports departures from the named preset.

| Setting | Initial value | Basis |
|---|---|---|
| DNA-binding arm | 36-42 nt | Ellenberg/LoopTrace Primary FISH probe design; also the main early-pilot range |
| Arm Tm | Strictly 42 < Tm < 50 C | Published numerical window; strict endpoint comparison from the existing implementation |
| Hybridization context | 2xSSC / 50% formamide | Reported in paper |
| Detailed Tm calculation | NN3, 390 mM Na, 2/0 nM, entropy salt correction 5, 0.65 C per percent FA | Explicit retained assumptions, not recovered exact paper settings |
| GC minimum / maximum | No independent minimum / 76% maximum | Retained project cap from documented G5 audit; not a verified earliest-pilot cutoff |
| Oligo spacing | At least 5 bp between extended footprints | User request retained; paper used 2 nt |
| Structure, off-target and RNA thresholds | Visible adjustable heuristics | Our QC rules, not a claim of the paper's LDA classifier |
| Layout | User-selected equal adjacent bins | The paper's 10 kb/up-to-100-primary example is separately selectable |

The paper's 40-50% GC restriction belongs to the short imaging oligos, not the primary
DNA-binding arms. The first project pool used mostly 36-42 nt candidates; four long
R2 probes were special additions, not permission to call 36-90 nt the paper default.
The broader G5-length option remains explicit. Neither preset changes old sequences.

The actual working references are GRCh38.p14, whereas the paper reports GRCh38.p12.
The program retains the user's simplified docking-TT-arm-TT-docking architecture,
not the paper's amplified multibarcode architecture. LDA and genome-wide k-mer
filtering are not implemented in this program; that status remains visible.

## The strand rule, without confusing DNA and RNA

TOP (+) is always the supplied forward DNA sequence; BOTTOM (-) is its complement.
For a plus-sense gene, the RNA has the same sequence polarity as TOP. A primary
with the same sequence as TOP is not complementary to that corresponding RNA, but
it is complementary to BOTTOM DNA. For a minus-sense gene the situation reverses.
Only the binding-arm orientation changes; docking sequences remain complementary
to the selected imager and are not reverse-complemented together with the full oligo.

| RNA sense within a candidate interval | Automatic bound DNA strand |
|---|---|
| Plus only | BOTTOM (-) |
| Minus only | TOP (+) |
| Both senses or unknown annotation strand | No candidate in the conflicting span |
| No annotated/declared span | BOTTOM convention plus independent RNA screen; not proof of no transcription |

This rule acts at each candidate location and on EVERY exact repeat copy. It is not
just a visual label and not a post-selection warning. A forced opposite strand is
rejected when it conflicts with a gene rather than overriding the RNA safeguard.

## What is required from a construct user

Each synthetic/local reference begins as UNKNOWN. The UI asks for forward, reverse,
both, reviewed annotation, manual transcript intervals, or no known transcription.
A separate checkbox confirms the review. A single declared direction conservatively
covers the entire sequence; bounded manual spans allow multiple genes and directions.
Known annotations are retained, not silently erased by a user declaration. Unknown
feature directions exclude both orientations over the feature span.

The declaration concerns transcription relative to the supplied sequence, not genomic
plus/minus guessed from a gene name. Include introns and untranslated sequence in
transcript spans. A coding sequence alone is not a full RNA model. Opposite-strand
overlaps may reduce the number of designable primaries, and shortages are reported.

## Independent RNA-complementarity checks remain

The direction guard is combined with the existing RNA search, strengthened in 1.1:
local unspliced gene spans are generated as RNA-sense references; spliced annotated
mRNA/other available RNA features are also searched. The complete primary and arm
are screened against supplied sequences and the installed human RNA reference in
full mode. The default excludes exact complementary 18-base tracts or the configured
long/high-identity RNA near matches. These are sequence heuristics, not RNA-DNA Tm.

Human DNA on-target exemptions and the human-DNA report-only override do not excuse
RNA matches. RNA avoidance requires RNA screening to stay on. Disabling avoidance
is explicit and prominently recorded. Unknown gene/transcript annotations and absent
RNA databases do not count as passes; full mode stops when required data are missing.

No finite reference set can establish absence of all antisense, enhancer, novel,
unspliced or read-through RNAs. Externally annotated public genes are considered
regardless of the selected cell line. H3K27ac remains a context track, not an RNA
absence test. Complete custom construct/transcript input is still needed for its
own comprehensive scope. Aqueous folding and heuristic BLAST limits are unchanged.

## Viewing and reproducibility

The dark interface covers the editor, QC controls, charts, sequence inspector,
source metadata and glossary. New results save the transcription declarations and
parameter provenance, with per-primary strand evidence. Mixed-strand region counts
are reported without assuming both strands are exposed on the same RASER-treated
DNA copy. Historical inputs and order files remain byte-identical.

The bundled tutorial is synthetic generated DNA and a local-only selection example.
It is not the CFTR target and must not be substituted for the historical order.
