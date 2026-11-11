Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

# Sources and provenance for version 1.1

## User-provided primary sources

- Beckwith et al. Nanoscale 3D DNA tracing in non-denatured cells resolves the
  Cohesin-dependent loop architecture of the genome in situ. Nature Communications
  16, 6673 (2025). DOI: 10.1038/s41467-025-61689-y. Supplied paper PDF,
  Primary FISH probe design (page 11): 36-42 nt, 42-50 C accounting for 2xSSC/50%
  formamide, minimum 2 nt spacing, 42 C LDA model/0.8, 5 or 10 off-target 18-mers,
  up to 100 primaries in 10 kb. Imager design in the next paragraph has 40-50% GC;
  that is not a primary-arm criterion. Exact complete Tm command is not supplied.
- RASERFISH_primary_probe_design_parameters_and_logic.pdf: initial 68-oligo pilot,
  mostly short arms, four long R2 support exceptions, approximate GC/length Tm and
  simplified two-handle architecture. Does not supply a verified primary-GC gate
  for the earliest pool. Observed GC values must not be relabeled as selection rules.
- REVIEW_README.txt / later G5 audit: GC <=76%, no independent minimum, no homopolymer
  >=6, original extended length search. These are project fallback settings.
- Imaging-oligos-I-currently-have.txt: the actual six imager sequences and dye labels.
- Previous Design Studio v1.0 source ZIP, including existing_design/ and inputs/:
  the software basis for this update. Historical sequences are not altered.

## Primary technical documentation checked

- LoopTrace paper: https://www.nature.com/articles/s41467-025-61689-y
- Ensembl overlap API (supports gene, transcript and exon features and strand):
  https://rest.ensembl.org/documentation/info/overlap_region
- OligoMiner blockParse source (parameterized Tm, not proof of paper settings):
  https://github.com/beliveau-lab/OligoMiner/blob/master/blockParse.py
- Biopython thermodynamic API:
  https://biopython.org/docs/1.86/api/Bio.SeqUtils.MeltingTemp.html
- ViennaRNA DNA parameter sets and cofold interpretation:
  https://viennarna.readthedocs.io/en/latest/eval/energy_parameters.html
  https://viennarna.readthedocs.io/en/latest/man/RNAcofold.html
- UCSC API: https://genome.ucsc.edu/goldenPath/help/api.html
- ENCODE API: https://www.encodeproject.org/help/rest-api/

## Explicit project choices, not published validation

Gene-direction exclusions, the 18-base RNA-complementarity gate, whole-gene span
screening, five-base extended footprints, full-primary architecture, folding/pair
thresholds, balanced greedy selection and candidate sampling are custom documented
choices. The native-target allowlist is not evidence that human RNA is safe. None
of these sources experimentally validates the newly generated candidate sequences.


## Installation sources rechecked 14 October 2026
- Official Miniforge release 26.7.2-0 and platform-specific SHA-256 sidecar.
  https://github.com/conda-forge/miniforge/releases/tag/26.7.2-0
- ViennaRNA 2.7.2 PyPI wheels (Python 3.11 macOS ARM64/x86_64 supported).
  https://pypi.org/project/ViennaRNA/2.7.2/
- NCBI Datasets accession download, --include genome,rna,seq-report.
  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/
The installation happens on the user's computer; these source checks do not
mean a macOS installation or large genome download was executed in the build.
