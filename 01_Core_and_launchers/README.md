# RASERFISH Design Studio

Local browser application for designing and reviewing RASER-FISH primary probes.
Version **1.6.1**, source snapshot **11 November 2026**, taken from the current
installed v1.6 folder. Includes the static chromatin-tracing background, probe
coverage views, bounded BLAST processing, compact outputs and checkpoint recovery.

## Start with an existing installation

Keep the existing environment and databases at their original location. From this
source folder, point the launcher to that installation:

```sh
bash start_studio.sh --qc-root /absolute/path/to/existing/Studio
```

If this source is already in the installed folder with `.qc_env` and `references`, use:

```sh
bash start_studio.sh
```

Open the address printed in Terminal, normally `http://127.0.0.1:8769/`. Keep
Terminal open. The launcher reuses an existing matching Studio and uses another
free local port when needed. Launching does not submit a new design; previously
enrolled, eligible interrupted work may resume from verified checkpoints.

## First installation on another machine

All four upload chunks must be combined into one repository, with `studio.py`,
`web/`, `engine/`, `inputs/`, `setup/` and `tests/` at the root.
The supplied installer supports macOS and Linux; this snapshot was checked on
macOS using an existing Python 3.11 environment, not a fresh installation.

Keep the new installation at its final location, with no spaces in the path.
The following command installs dependencies and downloads/indexes public human
genome and annotated-RNA databases. Plan roughly 30 GB free, with more needed for
working data. It is unnecessary for an already installed copy.

```sh
bash install.sh --no-open
bash start_studio.sh
```

Python packages alone are insufficient: BLAST+, its reference indexes and the
ViennaRNA DNA folding engine are also required. JavaScript tests additionally
use Node on PATH. No installation or reference download was performed while
creating this source package.

## Design policy

The ERRFI1 default includes a compact, public GRCh38 sequence and annotation
context: chr1:7,820,001–8,220,000, with the default DP1 design window at
chr1:8,009,019–8,029,018 on Top (+).

- Strict 45 < Tm < 55 °C; OligoMiner-style 2 nM profile; GC 0–76%.
- Preferred arms 36–42 nt, maximum 90 nt; 5-bp extended-footprint spacing.
- Genomic screening and annotated-RNA screening enabled; whole-gene exclusion off.
- High genomic off-target cutoff 42 °C, allowance 0; moderate cutoff 35 °C,
  allowance 10; unresolved long hits rejected; strong RNA-risk allowance 0.
- Individual folding/accessibility and self-dimer checks retained in full mode;
  minimum arm accessibility 0.40; folding at 42 °C and 0.39 M sodium.
- Cross-primary interaction screening defaults **OFF** and is reported as
  **NOT RUN - USER DISABLED**, within and between regions. Orders remain **DRAFT**.
- Spatially balanced search uses 400 preferred + 200 extended candidates per
  region. Maximize sites; benchmark 100 and desired minimum 50.
- One active design, two BLAST threads, bounded processing, checkpoint identity
  checks and memory/disk guards. A stopped analysis remains incomplete.

The November optimization groups four candidates per external BLAST process while
bounding native input groups. It does not remove queries, hits or QC thresholds.
The code uses CPU engines; this snapshot does not add a GPU backend.

These are project heuristics, not a reproduction of the complete published
OligoMiner pipeline or experimental validation. Annotation does not establish
expression; aqueous folding does not establish in-cell behavior. Review saved
scope and QC before using any output. See `METHODS_AND_LIMITATIONS.md`, the dated
release notes and `PERFORMANCE_UPDATE_20261103.txt` for provenance and limitations.

## What is included

Application source, UI, engine adapters, setup scripts, default inputs, synthetic
tutorial fixtures, regression tests and selected scientific baseline records.
The static background is embedded in `web/styles.css`, so it needs no network
request or animation. Its artwork comes from the project owner's website,
https://chb-bohrer.github.io/.

The environment, full search databases, real experimental runs, oligo orders,
saved personal configurations, caches, logs and backups are excluded. `.gitignore`
keeps these local files out of future commits. The small files in `inputs/` and
the synthetic examples are required source fixtures and are intentionally included.

## Checks

From the repository root, with the required tools on PATH:

```sh
PYTHONDONTWRITEBYTECODE=1 .qc_env/bin/python -m unittest discover -s tests -v
node --test tests/*.js
```

Alternatively use the Python executable from an existing installation. The
optional parity check against private six-target results skips when those files
are absent. See `PACKAGE_VALIDATION.txt` for the checks actually run on this package.

`PACKAGE_SHA256.json` records the packaged files, source hashes and the few
packaging-only edits. Runtime Python, browser assets, engine code and scientific
defaults were copied unchanged. No software license has been added by this
packaging step; dependencies retain their respective licenses.
