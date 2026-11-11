Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

# RASERFISH Design Studio 1.6.1

Updated again 16 October 2026 with verified checkpoints and bounded automatic recovery; see `AUTOMATIC_RECOVERY_2026-10-16.txt` for current behavior and validation. The earlier BLAST reliability repair remains in place; see `REPAIR_REPORT_2026-10-16.txt` for current validation. The original 15 October validation below describes the preceding build.

15 October 2026

This is a source update to the verified installed `RASERFISH_Design_Studio_v1_6` folder. It preserves the existing environment, reference databases, saved settings and design runs. Changed source files are backed up before replacement. No reinstall or database download is required.

## Provisional interaction policy

**Probe-probe interaction screening defaults to OFF.** Cross-primary exact-complementarity/stem searches and thermodynamic pairwise cofolds are skipped during selection, spatial repair and final export. This applies both within regions and between regions; it is not a raised threshold or an analysis performed and then hidden.

Individual-primary folding/accessibility and self-dimer checks remain active in full mode. Self-dimer records are explicitly separate from cross-primary records. Reports and exports mark omitted cross-primary analysis **NOT RUN - USER DISABLED**, never PASS. Order files remain **DRAFT**. Turning screening ON requests the original cross-primary checks and thresholds. Historical saved results are interpreted using their saved scope; missing old result fields are not retroactively labeled disabled.

## Preserved ERRFI1 settings

Top (+) throughout; genomic and annotated-RNA screening ON; whole-gene RNA exclusion OFF; strict **45 < Tm < 55 °C** with the OligoMiner-style 2 nM profile; GC 0–76%; adaptive arms 36–42 nt preferred and 90 nt maximum; 5-bp extended-footprint spacing. High off-target cutoff 42 °C with allowance 0; moderate cutoff 35 °C with allowance 10; unresolved long genomic hits rejected; strong RNA-risk allowance 0. Minimum arm accessibility 0.40; DNA folding at 42 °C and 0.39 M sodium.

Self-dimer rejection remains exact stem **>11 bp**, or cofold gain **<−12 kcal/mol** when intermolecular pairs exist; equality remains allowed. Maximize-sites selection retains benchmark 100 and desired minimum 50, the spatially balanced 400 preferred + 200 extended shortlist, spatial repair, and compact outputs. These scientific values were not relaxed.

## Memory and execution

- BLAST defaults to **two threads** and is capped at two. Importing an older higher thread value reduces it to two. A new configuration without the interaction setting defaults OFF; an explicit saved ON setting remains ON.
- A shared local lock permits one active design across updated app instances for the same macOS user, independent of folder or port.
- BLAST parsing and specificity aggregation use a disk-backed alignment store. Retained evidence is not truncated by an alignment-count cap. Compact exports omit large alignment collections before traversing/copying results. Full audit can still require substantial disk space.
- Folding caches are bounded. Native fold compounds and base-pair probability matrices are released after each calculation. Individual folding and genomic BLAST remain guarded even with cross-primary checks disabled.
- At the user's request, **available-memory reserve enforcement is disabled**, both at startup and during a design. The former 4 GiB startup and 2.4 GiB running reserve cutoffs no longer block or stop a run. Estimated available memory remains in the logs, with `available_memory_enforcement: disabled_by_user`. On this **16 GiB Mac**, the guard limits the design process plus children to **8 GiB resident memory or 8 GiB macOS physical footprint**, including charged compressed/swapped memory. The footprint check was added on 16 October after macOS reports demonstrated memory missed by RSS. The available-memory estimate uses free pages plus the smaller of inactive and file-backed pages. The supervisor samples approximately every five seconds, including during native work. Missing memory metrics block or stop analysis.
- The 16 October repair runs **one candidate / two arm-and-full-primary queries per BLAST process** and parses each successful batch into disk-backed evidence before launching the next. The **600-candidate budget is unchanged**, as are all BLAST sensitivity flags, hit limits and intended-site rules. Only successful new temporary outputs are removed in Standard mode; interrupted outputs are preserved and never reused. Every source remains incomplete until all its batches complete.
- Subprocess and supervisor journals preserve executable, arguments, working directory, PID/group, actual return code or signal, separate stdout/stderr, and explicit stop reasons. Owned workers are stopped on resource failures or lost supervision. Large downloads are streamed. Failed exports remain incomplete.
- New run folders retain `memory_usage.jsonl`, `worker_memory_usage.jsonl` and lifecycle status for diagnosis. A stopped, failed, cancelled or unfinished run remains **INCOMPLETE**; partial files are not served as completed results or orders. A completed provisional run means its requested analyses finished, not that omitted analyses passed.

The guard uses estimates and periodic observations. It reduces risk but cannot guarantee protection from a rapid allocation or operating-system termination. Disabling cross-primary checks alone is not treated as the memory fix.

## Launch

From the existing installed folder, open `Open_Studio.command` or run:

```sh
bash "<local-path-omitted>"
```

The launcher reuses the existing environment and databases, reuses an already open Studio, and otherwise prefers port 8769 with a free-port fallback. Use the actual address printed in Terminal. It restarts unexpected server exits at most twice. New Standard runs can automatically resume eligible failures from verified checkpoints, at most twice after 30 and 120 seconds. Launching does not start a new design; an explicitly enrolled interrupted design may resume. Cancel, intentional shutdown and safety/evidence failures stop automatic recovery.

Validation results for this update are recorded separately. Historical v1.6 test counts and example runs do not establish a new completed ERRFI1 design or experimental validation.

## Verified checks

340 Python regression tests and 16 JavaScript pairing tests passed. The tiny full-engine test used real genome/RNA BLAST fixtures and ViennaRNA, proving zero cross-primary calls with retained monomer and self-dimer checks. Two real-folding fixtures matched the actual installed v1.6.0 outputs exactly. Process tests exercised cross-instance locking, cancellation, simulated memory stops and incomplete exports. No full design was started. Details are in `verification/release_validation_v161.json`.
