Packaging note (2026-11-11): This is a dated development record. See README.md for current setup and PACKAGE_NOTES.md for this source snapshot. Paths to private evidence have been omitted; those files are not included.

> Historical feature notes. For the CURRENT 1.5 defaults, installation and RNA/strand policy, see START_HERE.txt and METHODS_AND_LIMITATIONS.md. Strict RNA exclusion is no longer the default.

# Design Studio 1.4: live imager / docking base-pairing guide

This is a display-only update to Design Studio 1.3. The scientific candidate
selection, Tm/GC, RNA-direction, off-target, folding, pair and spacing code and
settings are unchanged. No experimental probe set was redesigned or reordered.

## Where to see it

- Open `OPEN_BINDING_GUIDE.html` to go directly to Step 3 in a standalone preview.
- In the live app, Step 3 has a **Follow the paired bases** panel. Choose a readout
  in its preview dropdown or click **Show binding** in a region/inventory row.
  Changing the actual region assignment also updates the guide.
- Adding a new imager OR its docking site displays the resolved duplex before
  confirmation. Invalid input clears the earlier drawing; the diagram does not
  bypass the existing confirmation requirement or imply a biological QC pass.
- Step 6's candidate inspector shows the SAVED candidate's pair and architecture.
  Changing an unsaved dropdown never changes that saved diagram or sequence.
- New standalone design-result HTML includes the rendering code and works offline.
  Older exported HTML files are not altered by extracting this version.

## What the diagram means

Both stored sequences are always written 5-prime to 3-prime. The purple docking
row is drawn 5-prime to 3-prime from left to right. The gold imager row is drawn
3-prime to 5-prime so the two complementary sequences face one another.
Only the DISPLAY order of the imager is reversed; no ordered sequence is changed.

Recorded Dp7 example:

    Imager in paired display:  3' ATAAGCATAGCC 5'
                                 ||||||||||||
    Primary docking site:     5' TATTCGTATCGG 3'

    Imager for copying/order:  5' CCGATACGAATA 3'
    Docking for primary:      5' TATTCGTATCGG 3'

Copy sequences from the explicitly labeled canonical 5-prime-to-3-prime fields,
not the reversed gold display row. Hover or keyboard-focus a base-pair column to
see its position in each stored sequence. Docking position 1 faces imager position
N. A pairs with T; C pairs with G.

The overview shows one continuous primary:

    docking site 1 - TT - DNA-binding arm - TT - docking site 2

Each of its two identical docking sites is shown with a SEPARATE imaging molecule.
The teal middle arm has a different job: recognizing the genomic/construct target.
The schematic is not to scale. The dye label is shown, but its attachment end is
not inferred from the inventory. Two occupied handles illustrate potential pairing,
not a measured occupancy, folding result, imager Tm, fluorescence or QC pass.
This view does not show or redefine the top/bottom genomic DNA strand.

## Data integrity and implementation

`web/imager_binding.js` generates the SVG directly from the actual sequence pair.
There is no static decorative base sequence, external drawing service or image
library. Every valid input length (8-40 nt) is drawn without truncating bases.
Unsupported/inconsistent pairs return an explicit error, not a plausible drawing.
A saved candidate's full primary must agree with docking-TT-arm-TT-docking before
its assembly diagram is drawn. Labels are HTML/SVG escaped.

Pairing SVG and primary-schematic SVG can be exported separately. Exported SVGs
include their own colors, strand labels and descriptive text. The pairing SVG also
records the canonical imager and docking sequences in its attributes/description.
The buttons export pictures only, not a probe order. No external network access is
required for rendering or SVG export.

The preview dropdown is DISPLAY ONLY and does not assign a readout to a region.
Use the actual region-assignment dropdown to change an assignment. Such a sequence
change retains the existing requirement to rerun design/QC before a new order.

## Testing

See TESTING.txt and testing/BINDING_UI_RESULTS.json. New JavaScript tests cover
all six supplied imagers, all accepted lengths, antiparallel positions, inconsistent
pairs, user-text escaping, primary-assembly consistency and export metadata.
Chromium checks cover the live UI, offline preview, new-entry form, invalid edits,
real dropdown changes, mouse/keyboard interaction, SVG downloads, narrow viewports,
and saved-result export/reopening. Scientific/input file hashes are compared with
the supplied 1.3 archive in testing/INPUT_PRESERVATION_1_4.json.

No genome searches or full-folding analyses were needed or rerun for this UI change.
