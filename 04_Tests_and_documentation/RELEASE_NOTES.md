Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

> Historical v1.6.0 release notes. For the installed v1.6.1 provisional policy and memory safeguards, see [RELEASE_NOTES_1_6_1.md](RELEASE_NOTES_1_6_1.md). Statements below about all-pair screening, four threads and port 8769 describe v1.6.0.

# RASERFISH Design Studio 1.6.0

15 October 2026

This release makes the verified ERRFI1 / DP1 run the startup configuration, distributes the full-QC shortlist across the target, and reduces retained output size. Scientific acceptance thresholds and their exact comparison operators remain unchanged from the saved ERRFI1 run. It is a complete app folder; do not merge it into an older installation.

## Start with the existing installation

Extract `RASERFISH_Design_Studio_v1_6.zip` beside your installed `RASERFISH_Design_Studio_v1_5` folder. Open `Open_Studio.command`, or run `bash start_studio.sh` inside the new folder. The launcher reuses the sibling v1.5 software and reference databases, starts v1.6 on port 8769, and saves new results under v1.6's `design_runs/`. If the existing installation is elsewhere, use `bash start_studio.sh --qc-root /path/to/installed/Studio`.

No new download or reinstall is required when reusing the working installation. Fresh installations can still use `install.sh`. Keep an existing run running in its original folder; avoid starting another resource-heavy full run until it finishes. The old app, project sources and run outputs were not changed during this build.

## Verified ERRFI1 defaults

The reference was recovered from actual run `20261015_100112_62f902`, with its saved config, parameter provenance, generation counts and selected-site tables. The matching numbered settings downloads corroborate the settings. An older unsuffixed settings file had different values and was not used as the default source.

| Setting | v1.6 default |
|---|---|
| Reference context | GRCh38 chr1:7,820,001–8,220,000; 400,000 bp |
| Target | chr1:8,009,019–8,029,018; local 189019–209018; one 20,000 bp region; strict tail policy |
| Imager | DP1, 5′-CTAACGAATCGG-3′; docking 5′-CCGATTCGTTAG-3′ |
| Mode / objective | Full QC / maximize compatible exact binding sites; target benchmark 100, minimum benchmark 50 |
| DNA-bound strand | Top (+), trace-wide |
| On-target Tm | Strictly 45 < Tm < 55 °C; OligoMiner-style 2 nM profile |
| Arm lengths | Adaptive, 36–42 nt preferred; extensions through 90 nt |
| GC / homopolymer | 0–76%; maximum homopolymer 5 |
| Off-target cutoffs | High 42 °C, allowed 0; moderate 35 °C, allowed 10; unresolved long hits rejected |
| RNA | Annotated human and supplied transcript screening ON; strong-risk allowance 0; exact-tract setting 18; whole-gene exclusion OFF |
| Folding | 42 °C; sodium 0.39 M; minimum arm-unpaired probability 0.4 |
| Pair gates | Reject cofold gain below −12 kcal/mol when intermolecular pairs exist, or exact complementary stem above 11 bases; equality at either limit remains allowed |
| Footprints | Extended; minimum gap 5 bp; repeated exact sites allowed only within their assigned region |
| Human hits | Intended intervals exempted by coordinates and strand; unresolved bins not allowed |
| Search controls | Step 1 bp; 2 variants/start; 5,000 total candidate cap; 4 threads |
| Output | Standard (compact) |

`inputs/errfi1_defaults.json` is the shared source for the active ERRFI1 preset. Startup, field defaults, preset application and full ERRFI1 restore agree. Applying the policy preset restores its verified thresholds and RNA/strand choices. “Restore ERRFI1 setup” also restores target bounds and DP1. Prior temperature-comparison presets are retained as historical code/documentation, outside the active preset menu. An explicitly imported custom configuration retains its requested settings.

## Spatial candidate search

The only changed default search budgets are **400 preferred-length + 200 extended candidates per region**, replacing 160 + 80. The target is divided into fixed 1 kb internal sectors; candidates are selected across sectors and dispersed within each sector after the same local filters. This gives 20 internal sectors for the 20 kb ERRFI1 window; it does not create 20 separate imager regions.

After full QC, selection prioritizes compatible exact-site count, then coverage gaps and distribution. A bounded deterministic repair step considers alternate passing candidates and can replace up to two blocking selections before refilling. It checks the same footprint, strand, individual and pair gates; it never accepts a repair that reduces selected-site count. Selection remains a finite heuristic, not a proof of the globally largest set.

Reports include 1 kb sector counts, largest uncovered gap including region edges, and median consecutive site-center spacing. The existing ten-bucket coverage statistic is retained for compatibility. The 100/50 benchmarks are reporting targets, not automatic QC relaxation or a cap in maximize mode.

On the actual ERRFI1 reference, a fresh candidate-generation check produced **600 unique Top (+) candidates: exactly 30 in each of 20 sectors** (20 preferred + 10 extended). Every pre-shortlist rejection/filter count matched the previous run. These candidates have not undergone a new full ERRFI1 BLAST/RNA/folding run. **This release does not claim more than the previous 41 selected probes.** The saved run has selected-site and pair tables; its current directory lacks a final result JSON/viewer/completion record, so the old run is not represented as a fully packaged completed result.

## Smaller saved results

**Standard**, the default, saves selected probes/sites and QC summaries, compact result JSON, a standalone `viewer.html`, per-sector coverage, pair summary, flagged pair details, selected QC flags and compressed candidate decisions. The viewer contains selected candidates only. Counts and rejection summaries still describe all reviewed candidates.

`run_config.json` is now a versioned manifest with filenames, lengths and SHA-256 hashes. Compressed FASTA/annotation assets store the actual inputs outside JSON. The 400 kb context is retained once as a compressed input asset when needed for reproducible local exact-site screening; it is absent from result JSON and the result viewer. The target slice is also saved explicitly. Settings export is a portable ZIP of this manifest and its assets; import checks the assets and restores the original inputs. A manifest alone is not a complete settings package. Legacy full-sequence settings JSON remains importable.

Raw BLAST output is still created temporarily while each database runs. After parsing and hashing, Standard deletes it from the new run. **Full audit** instead compresses raw searches, candidate/alignment details and all selected-pair records under `debug/`. Full audit can still be large. No raw alignments are embedded in the default JSON or HTML, and the former multi-gigabyte review CSV is no longer generated before compression. Working memory and temporary disk use can still be substantial, especially with 600 candidates.

Near-limit reporting flags do not change acceptance: selected-pair rows are flagged within 1 base of the stem setting or within 1 kcal/mol of the cofold setting when intermolecular pairing is present. All selected pairs are still evaluated. Full audit retains all pair rows; Standard retains their totals and flagged rows. The browser retries temporary job-status fetch failures while the server continues working.

Format regression on the existing ERRFI1 inputs reduced the 873,046-byte config to a 38,728-byte manifest and **171,421 bytes for the entire manifest/input bundle** (80.4% smaller combined). Candidate decisions compressed from 98,183 to 13,563 bytes. The 861 saved selected-pair rows yielded 7 reporting flags and zero failures. These are format-conversion measurements, not a new scientific rerun. Prior large raw files were measured by filesystem size only and left untouched.

## Verification and limits

The full regression suite passed **294 tests**, including real BLAST/ViennaRNA adapter tests, scientific source-hash checks, exact threshold boundaries, spatial balancing/repair, compact outputs, and settings round trips. The JavaScript pairing suite passed all **16 tests**. The actual ERRFI1 generation check and synthetic local-only end-to-end design also passed. The synthetic demo selected 46 probes across three regions; it is clearly labeled local-only and is not an ERRFI1 result.

`verification/` and `validation/` contain the evidence, with source hashes for the unchanged scientific routines/modules. `TESTING.txt` records scope. A new full ERRFI1 QC run, fresh software/reference installation and experimental performance validation were not performed. The five-window plan remains planning material and still needs separate runs plus combined-pool QC. No acceptance threshold, reference database or thermodynamic model was silently relaxed.
