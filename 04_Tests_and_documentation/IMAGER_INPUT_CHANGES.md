Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

> Historical feature notes. For the CURRENT 1.5 defaults, installation and RNA/strand policy, see START_HERE.txt and METHODS_AND_LIMITATIONS.md. Strict RNA exclusion is no longer the default.

> Historical 1.3 input-convention notes. Version 1.4 retains the same stored
> sequences but replaces the old text alignment with a visual diagram whose
> IMAGER row is reversed for display (3-prime to 5-prime) above the canonical
> docking row (5-prime to 3-prime). See BINDING_DIAGRAM_CHANGES.md.

# Design Studio 1.3: imaging oligo OR its docking sequence

## User-facing behavior

Step 3 has a required choice: **Imaging oligo itself** or **Sequence it binds to:
docking site on the primary**. Both are entered in their own 5-prime-to-3-prime
orientation. There is no default polarity guess for a new entry.

The program previews both sequences before saving, labels which is the separate
fluorescent oligo and which is appended to the primary, and requires a review
checkbox. Editing the sequence/type clears the confirmation. An invalid edit
clears the previous preview instead of leaving an apparently valid stale pair.

The optional antiparallel diagram shows the docking strand 3-prime to 5-prime
ONLY for the base-pairing display. The stored and entered docking sequence is
still 5-prime to 3-prime. Direction labels and dye notation are not sequence input.

For the recorded Dp7 pair:

```
Imaging oligo 5'-CCGATACGAATA-3'
Docking site  5'-TATTCGTATCGG-3'

Primary 5'-TATTCGTATCGG-TT-[DNA-binding arm]-TT-TATTCGTATCGG-3'
```

The docking site is NOT the genomic/construct target sequence. The latter is
loaded in step 1. Changing entry convention for an equivalent sequence pair
must not change a primary or any QC value. Changing the physical pair requires
new full-primary checks, as before. This input validation does not establish
barcode orthogonality or experimental imager performance.

## Internal representation and API

Legacy input is unchanged:

```json
{"id":"MY_DP7","sequence":"CCGATACGAATA"}
```

New explicit entries may use either of these forms:

```json
{"id":"MY_DP7","input_kind":"imager","input_sequence_5to3":"CCGATACGAATA"}
```

```json
{"id":"MY_DP7","input_kind":"docking","input_sequence_5to3":"TATTCGTATCGG"}
```

Both resolve to `sequence = CCGATACGAATA` and `docking = TATTCGTATCGG`.
Explicit alias fields are `imager_sequence_5to3` and `docking_sequence_5to3`.
The original normalized entry is retained as input_kind/input_sequence_5to3.
The generic `sequence` field ALWAYS means the imaging oligo, never the raw
docking entry. If redundant fields disagree, validation stops. New ambiguous
input_sequence_5to3 without a type also stops. Stored orientation is `5to3`.

The browser and `/api/imagers/resolve` validate the same conventions. The
`/api/design` endpoint validates before starting a job, and design() validates
again. Uppercase and whitespace normalization is permitted; non-ACGT letters,
numbers, dye notation, or direction labels are rejected. Existing 8-40-nt input
limits are unchanged. This is an input bound, not a recommended imager length.
Duplicate canonical imager/docking pairs are rejected across both modes.

## Outputs and persistence

- Step-3 inventory and assignment tables show both sequences.
- Candidate inspectors use the SAVED candidate pair, not unsaved dropdown edits.
- New selected_primaries.csv has Imager_Sequence_5to3,
  Docking_Sequence_in_Primary_5to3, Imager_Input_Type, and
  Imager_Entered_Sequence_5to3, separately from DNA_Binding_Arm_5to3 and
  Full_Primary_5to3.
- New imager_inventory.csv and readout_assignments.csv include entered convention
  and bases; they are reference tables, not primary-oligo orders.
- Settings/result JSON and offline results preserve custom pair provenance.
- Additions live in the current session and saved settings/results, not a global
  on-disk inventory that automatically persists between browser sessions.

The existing IDT order format stays Pool name / Sequence, with the same complete
primary sequence for equivalent readout pairs. No old order file is rewritten.

## Tests and limitations

178 automated tests passed; two pre-existing BLAST/ViennaRNA engine tests were
skipped because those engines are not installed in this execution environment.
106 browser/API/standalone checks passed (38 previous UI checks, 10 previous
preview checks and 58 new entry-convention checks). The new tests include an
actual local-only design with each entry type; equivalent entries produced
identical candidate primaries, selected IDs and draft ordering CSV bytes.

The interface was inspected in Chromium, including a narrow/mobile viewport.
External genome/RNA searches, live annotation retrieval and full folding were
NOT rerun for this input/reporting-only update. Their earlier modeling and
reference-coverage limitations remain unchanged. The existing synthetic tutorial
is still a local-only demonstration, not a new CFTR replacement order.
