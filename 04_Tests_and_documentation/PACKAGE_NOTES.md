# Source package — 11 November 2026

The source was copied from the installed `RASERFISH_Design_Studio_v1_6` folder.
Its actual version is 1.6.1. It was not reconstructed from an earlier release.

Runtime Python, engine and setup code, launch scripts, browser assets and default
input files are unchanged. The package includes the accepted static chromatin
illustrations, not the superseded animation or nucleus illustration.

Packaging-only edits:

- Added a current README, quick-start text, ignore rules, package notes and checksums.
- Removed personal absolute workstation paths from copied documentation and
  historical baseline metadata. No scientific values or sequences were changed.
- Made two test files portable: Node is discovered on PATH, and the optional
  private saved-run comparison looks relative to this source folder.
- Updated the old source-preservation test to recognize the November 3 cache
  extraction. It verifies the exact helper/wrapper identities, then normalizes
  only those substitutions and still compares against the original baseline
  hashes. Existing behavior-parity tests remain in force. The application code
  and original scientific baseline manifest were not altered.
- Marked dated release/development notes as historical. Their old test counts,
  memory observations and benchmarks describe those dates, not a fresh run today.
- Kept only the synthetic tutorial JSON used by the app and regression tests;
  omitted its duplicate exports and stale standalone HTML. These tutorial
  results retain their original saved analysis policy and are not new probe orders.
- Retained three small baseline verification records needed by defaults/provenance
  checks. Other machine-specific verification logs and snapshots were excluded.

Excluded directories include `.qc_env`, `.tools`, `references`, `design_runs`,
`experiments`, `existing_design`, `prepared_settings`, `plans`, `cache`, `logs`,
`setup_logs`, backups and old test-output folders. Historical real probe sets and
order files are not part of this software upload. Their installed originals remain
untouched. The old package checksum file was replaced with a new manifest for
this exact source snapshot.

All chunks must be uploaded before installing or launching a newly downloaded
repository. No automatic GitHub upload, publishing, installation, database download
or production design was performed.
