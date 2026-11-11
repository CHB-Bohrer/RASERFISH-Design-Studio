#!/usr/bin/env python3
"""Create a self-contained, offline primary-probe binding viewer.

Python 3.9+; standard library only. No genome download or web service is used.
All coordinates are 1-based, inclusive, relative to each supplied linear sequence.
The mapped sequence is the DNA-binding ARM, not the docking handles/spacers.

A = primary arm, S = forward reference:
    A == S[start:end]     -> physically binds the opposite (-) DNA strand.
    revcomp(A) == S[...]  -> physically binds the supplied (+) DNA strand.
Both searches are performed, retaining every full-arm exact occurrence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
from pathlib import Path
import re
import sys
from typing import Any, Iterable, Optional

HANDLES = {"DP7": "TATTCGTATCGG", "DP3": "CCGCGATATTAA"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
REFERENCE_ALPHABET = set("ACGTRYSWKMBDHVN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def dna(text: str, label: str, reference: bool = False) -> str:
    """Strip whitespace only; never silently remove annotation or invalid bases."""
    seq = re.sub(r"\s+", "", text).upper()
    alphabet = REFERENCE_ALPHABET if reference else set("ACGT")
    invalid = set(seq) - alphabet
    if not seq:
        raise ValueError(f"{label}: sequence is empty.")
    if invalid:
        raise ValueError(f"{label}: invalid character(s) {sorted(invalid)!r}. "
                         "Provide plain DNA sequence or FASTA, not annotated text.")
    return seq


def fasta_records(text: str, default_id: str = "R1") -> list[dict[str, str]]:
    """Read multi-FASTA or a single raw sequence. IUPAC letters are not wildcards."""
    text = text.lstrip("\ufeff")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("Reference file is empty.")
    if not lines[0].startswith(">"):
        return [{"id": default_id, "description": default_id,
                 "sequence": dna("".join(lines), default_id, reference=True)}]
    records: list[dict[str, str]] = []
    header: Optional[str] = None
    chunks: list[str] = []
    for line in lines + [">__END__"]:
        if line.startswith(">"):
            if header is not None:
                ident = header.split()[0]
                records.append({"id": ident, "description": header,
                                "sequence": dna("".join(chunks), ident, reference=True)})
            header = line[1:].strip()
            if not header:
                raise ValueError("FASTA record has an empty identifier.")
            chunks = []
        else:
            chunks.append(line)
    ids = [r["id"] for r in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Reference FASTA identifiers must be unique.")
    return records


def read_references(fasta: Optional[Path], references: Optional[list[str]],
                    input_dir: Path) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    if references:
        for spec in references:
            if "=" not in spec:
                raise ValueError("--reference must be ID=path, for example R1=gene.txt")
            ident, name = spec.split("=", 1)
            if not ident or any(c.isspace() for c in ident):
                raise ValueError("Reference IDs must be nonempty and contain no whitespace.")
            path = Path(name).expanduser()
            records = fasta_records(path.read_text(encoding="utf-8-sig"), ident)
            if len(records) != 1:
                raise ValueError(f"{path}: --reference requires one sequence per file; "
                                 "use --fasta for multi-FASTA.")
            rec = records[0]
            rec["id"] = ident
            rec["source_file"] = path.name
            result.append(rec)
    else:
        path = fasta or input_dir / "target_sequences.fasta"
        result = fasta_records(path.read_text(encoding="utf-8-sig"))
        for rec in result:
            rec["source_file"] = path.name
    if len({r["id"] for r in result}) != len(result):
        raise ValueError("Duplicate reference IDs.")
    for rec in result:
        rec["length"] = len(rec["sequence"])
        rec["sha256"] = hashlib.sha256(rec["sequence"].encode("ascii")).hexdigest()
    return result


def _field(row: dict[str, str], aliases: Iterable[str]) -> str:
    lower = {str(k).strip().lower(): (v or "").strip() for k, v in row.items() if k}
    for alias in aliases:
        value = lower.get(alias.lower(), "")
        if value:
            return value
    return ""


def parse_probe_csv(text: str, arms_only: bool = False) -> list[dict[str, Any]]:
    """Accept the supplied order CSV, the audit CSV, or an explicit arm table.

    An unrecognized Sequence column is NOT silently treated as an arm unless
    --arms-only is explicitly provided. Handles are removed only when both ends
    have the documented handle-TT / TT-handle architecture.
    """
    text = text.lstrip("\ufeff")
    if not text.strip():
        raise ValueError("Probe CSV is empty.")
    try:
        delimiter = csv.Sniffer().sniff(text[:8192], delimiters=",\t;").delimiter
    except csv.Error:
        delimiter = ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    if not reader.fieldnames:
        raise ValueError("Probe CSV has no header row.")
    probes: list[dict[str, Any]] = []
    names: set[str] = set()
    for row_number, row in enumerate(reader, 2):
        if not any((v or "").strip() for v in row.values() if isinstance(v, str)):
            continue
        if None in row:
            raise ValueError(f"Probe CSV row {row_number}: more values than column headers.")
        name = _field(row, ["Name", "Probe_Name", "probe_id", "id"])
        if not name:
            raise ValueError(f"Probe CSV row {row_number}: missing Name or Probe_Name.")
        if name in names:
            raise ValueError(f"Duplicate probe name: {name}")
        names.add(name)
        explicit_arm = _field(row, ["DNA_Binding_Arm_5to3", "Arm_5to3", "arm",
                                     "binding_arm", "DNA_binding_insert"])
        raw_sequence = _field(row, ["Full_Primary_5to3", "full_sequence", "Sequence", "full"])
        declared_readout = _field(row, ["Readout", "readout", "DP", "imager"]).upper()
        full = dna(raw_sequence, f"{name} full sequence") if raw_sequence else ""
        detected_readout = ""
        detected_arm = ""
        prefix = suffix = ""
        if full and not arms_only:
            for readout, handle in HANDLES.items():
                left, right = handle + "TT", "TT" + handle
                if len(full) > len(left) + len(right) and full.startswith(left) and full.endswith(right):
                    detected_readout = readout
                    detected_arm = full[len(left):-len(right)]
                    prefix, suffix = left, right
                    break
        if explicit_arm:
            arm = dna(explicit_arm, f"{name} arm")
            if detected_arm and arm != detected_arm:
                raise ValueError(f"{name}: explicit arm disagrees with the full oligo architecture.")
            if full and not detected_arm:
                pos = full.find(arm)
                if pos < 0:
                    raise ValueError(f"{name}: explicit arm is not contained in full oligo.")
                prefix, suffix = full[:pos], full[pos + len(arm):]
        elif arms_only and full:
            arm, full = full, ""
        elif detected_arm:
            arm = detected_arm
        else:
            raise ValueError(f"{name}: cannot safely identify the binding arm. "
                             "Supply DNA_Binding_Arm_5to3, use the documented DP7/DP3 "
                             "full architecture, or use --arms-only for an arm-only Sequence column.")
        if declared_readout and detected_readout and declared_readout != detected_readout:
            raise ValueError(f"{name}: declared readout disagrees with the docking handles.")
        readout = detected_readout or declared_readout or "Other"
        region = _field(row, ["Region", "expected_region", "Assigned_Region"])
        if not region:
            match = re.match(r"^(R\d+)(?:_|\b)", name)
            region = match.group(1) if match else ""
        probes.append({"name": name, "readout": readout,
                       "pool": _field(row, ["Pool"]), "expected_region": region,
                       "arm": arm, "full": full, "prefix": prefix, "suffix": suffix,
                       "arm_length": len(arm), "full_length": len(full) if full else None,
                       "gc_percent": 100.0 * (arm.count("G") + arm.count("C")) / len(arm)})
    if not probes:
        raise ValueError("Probe CSV contains no sequence rows.")
    return probes


def exact_starts(reference: str, needle: str) -> Iterable[int]:
    """Yield ALL zero-based starts, including overlapping/repetitive occurrences."""
    if not needle:
        raise ValueError("Cannot search for an empty arm.")
    start = reference.find(needle)
    while start >= 0:
        yield start
        start = reference.find(needle, start + 1)


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge 1-based inclusive intervals, counting overlapping bases only once."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def build_dataset(references: list[dict[str, Any]], probes: list[dict[str, Any]],
                  title: str = "CFTR DP7 / DP3 - minimum 5-base gaps") -> dict[str, Any]:
    hits: list[dict[str, Any]] = []
    warnings: list[str] = []
    for p_index, probe in enumerate(probes):
        probe["index"] = p_index
        probe["hits"] = []
        for ref in references:
            # A forward identity match implies binding to the OPPOSITE DNA strand.
            for bound_strand, needle in [("-", probe["arm"]), ("+", reverse_complement(probe["arm"]))]:
                for pos in exact_starts(ref["sequence"], needle):
                    hit = {"id": len(hits), "probe_index": p_index,
                           "reference_id": ref["id"], "start": pos + 1,
                           "end": pos + len(needle), "bound_strand": bound_strand,
                           "arm_match_to_forward": "identical" if bound_strand == "-" else "reverse_complement"}
                    hits.append(hit)
                    probe["hits"].append(hit["id"])
        if not probe["hits"]:
            warnings.append(f"No full-arm exact site for {probe['name']} in the supplied references.")
    unique_arms = {p["arm"] for p in probes}
    if len(unique_arms) != len(probes):
        warnings.append("Some CSV rows share the same binding arm. They are retained as separate named "
                        "probes, not deduplicated; site counts may therefore include duplicate oligo rows.")
    for ref in references:
        rh = [h for h in hits if h["reference_id"] == ref["id"]]
        ref["hits"] = [h["id"] for h in rh]
        ref["covered_intervals"] = merge_intervals((h["start"], h["end"]) for h in rh)
        coverage = sum(e - s + 1 for s, e in ref["covered_intervals"])
        gap_start = 1
        gaps = []
        for start, end in ref["covered_intervals"]:
            if start > gap_start:
                gaps.append({"start": gap_start, "end": start - 1, "length": start - gap_start})
            gap_start = end + 1
        if gap_start <= ref["length"]:
            gaps.append({"start": gap_start, "end": ref["length"], "length": ref["length"] - gap_start + 1})
        ref["gaps"] = gaps
        counts: dict[int, int] = {}
        for h in rh:
            counts[h["probe_index"]] = counts.get(h["probe_index"], 0) + 1
        ref["summary"] = {"unique_probe_rows": len(counts), "exact_sites": len(rh),
                          "covered_bases": coverage, "coverage_percent": 100.0 * coverage / ref["length"],
                          "repeat_probe_rows": sum(c > 1 for c in counts.values()),
                          "longest_uncovered_interval": max((g["length"] for g in gaps), default=0),
                          "opposite_strand_sites": sum(h["bound_strand"] == "-" for h in rh),
                          "forward_strand_sites": sum(h["bound_strand"] == "+" for h in rh)}
        if any(c not in "ACGT" for c in ref["sequence"]):
            warnings.append(f"{ref['id']}: ambiguous IUPAC bases are retained in coordinates but are "
                            "not wildcards; an exact A/C/G/T arm will not match across them.")
    return {"schema_version": 1, "title": title, "references": references, "probes": probes,
            "hits": hits, "warnings": warnings,
            "coordinate_system": "1-based inclusive; supplied linear sequence coordinates, not genomic coordinates",
            "mapping_method": "All full-arm exact matches in both DNA orientations; no partial or mismatched matching",
            "scope": "Sequence matching only. This is not experimental binding, accessibility, brightness, "
                     "RNA exclusion, or genome-wide off-target validation."}


def annotate_spacing(data: dict[str, Any], minimum_gap: int = 5) -> None:
    """Audit gaps between every adjacent full-arm site on each linear reference.

    Includes all repeat copies and both mapped strand orientations. This is
    conservative for mixed-strand inputs. It does not model partial binding or
    circular / unsupplied sequence junctions. Never hides a conflicting hit.
    """
    if minimum_gap < 0:
        raise ValueError("minimum_gap cannot be negative")
    data["required_min_gap_nt"] = minimum_gap
    data["spacing_audit"] = []
    for hit in data["hits"]:
        hit["gap_before_nt"] = None
        hit["gap_after_nt"] = None
        hit["previous_hit_id"] = None
        hit["next_hit_id"] = None
    for ref in data["references"]:
        hits = sorted((data["hits"][i] for i in ref["hits"]),
                      key=lambda h: (h["start"], h["end"], h["id"]))
        adjacent = []
        for a, b in zip(hits, hits[1:]):
            gap = b["start"] - a["end"] - 1
            a["gap_after_nt"] = b["gap_before_nt"] = gap
            a["next_hit_id"] = b["id"]
            b["previous_hit_id"] = a["id"]
            row = {"reference_id": ref["id"], "left_hit_id": a["id"],
                   "right_hit_id": b["id"], "start": a["end"] + 1,
                   "end": b["start"] - 1, "length": gap,
                   "required_min_gap_nt": minimum_gap, "pass": gap >= minimum_gap}
            adjacent.append(row)
            data["spacing_audit"].append(row)
        ref["adjacent_site_gaps"] = adjacent
        ref["summary"]["minimum_gap_nt"] = min((a["length"] for a in adjacent), default=None)
        ref["summary"]["spacing_pairs_checked"] = len(adjacent)
        ref["summary"]["spacing_violations"] = sum(not a["pass"] for a in adjacent)
        if ref["summary"]["spacing_violations"]:
            data["warnings"].append(f"{ref['id']}: {ref['summary']['spacing_violations']} adjacent "
                                    f"site pairs have fewer than {minimum_gap} intervening bases.")
    data["spacing_all_pass"] = all(row["pass"] for row in data["spacing_audit"])


def _write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(data: dict[str, Any], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    if "spacing_audit" not in data:
        annotate_spacing(data, 5)
    serialized = json.dumps(data, ensure_ascii=True, separators=(",", ":"))
    # Escaping '<' prevents an input-derived </script> from ending the data tag.
    safe_json = serialized.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    page = HTML_TEMPLATE.replace("__PAGE_TITLE__", html.escape(data["title"]))
    page = page.replace("__DATA_JSON__", safe_json)
    (out / "CFTR_probe_binding_viewer_DP7_DP3_gap5.html").write_text(page, encoding="utf-8")
    (out / "binding_map.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    fields = ["Probe_Name", "Readout", "Reference", "Start_1based", "End_1based_inclusive",
              "DNA_strand_bound", "Arm_match_to_forward_reference", "Arm_length_nt", "Arm_5to3", "Gap_before_nt", "Gap_after_nt"]
    rows = []
    for h in sorted(data["hits"], key=lambda h: (h["reference_id"], h["start"], h["end"], h["probe_index"], h["bound_strand"])):
        p = data["probes"][h["probe_index"]]
        rows.append(dict(zip(fields, [p["name"], p["readout"], h["reference_id"], h["start"], h["end"],
                                    h["bound_strand"], h["arm_match_to_forward"], p["arm_length"], p["arm"], h["gap_before_nt"], h["gap_after_nt"]])))
    _write_csv(out / "mapped_binding_sites_gap5.csv", fields, rows)
    spacing_fields = ["Reference", "Left_probe", "Left_start", "Left_end", "Right_probe",
                      "Right_start", "Right_end", "Gap_start", "Gap_end", "Gap_nt",
                      "Required_min_gap_nt", "Pass"]
    spacing_rows = []
    for item in data["spacing_audit"]:
        a, b = data["hits"][item["left_hit_id"]], data["hits"][item["right_hit_id"]]
        pa, pb = data["probes"][a["probe_index"]], data["probes"][b["probe_index"]]
        spacing_rows.append(dict(zip(spacing_fields, [item["reference_id"], pa["name"],
            a["start"], a["end"], pb["name"], b["start"], b["end"], item["start"], item["end"],
            item["length"], item["required_min_gap_nt"], "PASS" if item["pass"] else "FAIL"])))
    _write_csv(out / "spacing_audit_gap5.csv", spacing_fields, spacing_rows)
    _write_csv(out / "uncovered_intervals.csv", ["Reference", "Start_1based", "End_1based_inclusive", "Length_nt"],
               ({"Reference": r["id"], "Start_1based": g["start"], "End_1based_inclusive": g["end"], "Length_nt": g["length"]}
                for r in data["references"] for g in r["gaps"]))
    _write_csv(out / "mapping_summary.csv", ["Reference", "Length_nt", "Probe_rows", "Exact_sites",
                "Covered_bases", "Coverage_percent", "Repeat_probe_rows", "Opposite_strand_sites", "Forward_strand_sites", "Minimum_gap_nt", "Spacing_violations"],
               ({"Reference": r["id"], "Length_nt": r["length"], "Probe_rows": r["summary"]["unique_probe_rows"],
                 "Exact_sites": r["summary"]["exact_sites"], "Covered_bases": r["summary"]["covered_bases"],
                 "Coverage_percent": round(r["summary"]["coverage_percent"], 4),
                 "Repeat_probe_rows": r["summary"]["repeat_probe_rows"],
                 "Opposite_strand_sites": r["summary"]["opposite_strand_sites"],
                 "Forward_strand_sites": r["summary"]["forward_strand_sites"],
                 "Minimum_gap_nt": r["summary"]["minimum_gap_nt"],
                 "Spacing_violations": r["summary"]["spacing_violations"]} for r in data["references"]))


def main(argv: Optional[list[str]] = None) -> int:
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--probes", type=Path, default=base / "ORDER_ALL_DP7_DP3_gap5.csv",
                        help="Order CSV or explicit DNA-binding-arm CSV (defaults to bundled 163 gap-5 probes).")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--fasta", type=Path, help="Multi-FASTA with one record per reference sequence.")
    group.add_argument("--reference", action="append", metavar="ID=PATH",
                       help="Repeat for each individual raw DNA / single-record FASTA file.")
    parser.add_argument("--arms-only", action="store_true", help="Explicitly treat Sequence as an arm-only column.")
    parser.add_argument("--out", type=Path, default=base / "output", help="Output directory (default: output/).")
    parser.add_argument("--title", default="CFTR DP7 / DP3 - minimum 5-base gaps")
    parser.add_argument("--minimum-gap", type=int, default=5,
                        help="Required intervening bases between mapped sites (default: 5).")
    args = parser.parse_args(argv)
    try:
        refs = read_references(args.fasta, args.reference, base / "inputs")
        probes = parse_probe_csv(args.probes.read_text(encoding="utf-8-sig"), args.arms_only)
        data = build_dataset(refs, probes, args.title)
        if args.minimum_gap < 0:
            raise ValueError("--minimum-gap must be nonnegative.")
        annotate_spacing(data, args.minimum_gap)
        data["probe_source_file"] = args.probes.name
        data["probe_source_sha256"] = hashlib.sha256(args.probes.read_bytes()).hexdigest()
        write_outputs(data, args.out)
    except (OSError, ValueError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Mapped {len(probes)} probe rows to {len(data['hits'])} full-arm exact sites.")
    for ref in refs:
        summary = ref["summary"]
        print(f"  {ref['id']}: {ref['length']:,} nt, {summary['unique_probe_rows']} probes, "
              f"{summary['exact_sites']} sites, {summary['coverage_percent']:.2f}% sequence coverage")
    for warning in data["warnings"]:
        print(f"WARNING: {warning}", file=sys.stderr)
    print(f"Open in a browser: {(args.out / 'CFTR_probe_binding_viewer_DP7_DP3_gap5.html').resolve()}")
    return 0


# The complete offline viewer is embedded below; there is no CDN or dependency.
HTML_TEMPLATE = r"""
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light">
<title>__PAGE_TITLE__</title>
<style>
:root{--ink:#172b3b;--muted:#607180;--line:#dce3e8;--paper:#f4f6f8;--dp7:#176e7b;--dp3:#7654ac;--accent:#e9b550;--soft:#eef3f6}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.5 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
main{max-width:1540px;margin:auto;padding:30px 32px 36px}h1{font-size:32px;line-height:1.2;letter-spacing:-.7px;margin:7px 0 10px}h2{font-size:19px;margin:0}h3{font-size:15px;margin:16px 0 7px}p{margin:8px 0}.eyebrow{font-size:11px;letter-spacing:2px;font-weight:750;color:var(--dp7)}.lede{color:var(--muted);max-width:1020px}.headrow{display:flex;justify-content:space-between;gap:20px;align-items:flex-start}.badges{display:flex;gap:8px;flex-wrap:wrap;justify-content:flex-end;padding-top:5px}.badge{display:inline-block;background:white;border:1px solid var(--line);border-radius:20px;padding:5px 11px;font-size:12px;white-space:nowrap}.badge.good{background:#eaf3ef;color:#306749;border-color:#cddfd4}.note{font-size:13px;color:var(--muted)}.notice{background:#fff8e7;border-left:3px solid var(--accent);padding:10px 14px;border-radius:5px;margin:18px 0;font-size:13px}.overview{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin:20px 0}.card{background:white;border:1px solid var(--line);border-radius:12px;box-shadow:0 1px 2px #182a3910;min-width:0}.overview-card{padding:16px 18px 10px;cursor:default}.overview-card.active{border-color:#789da8;box-shadow:0 0 0 1px #789da830}.cardheading{display:flex;justify-content:space-between;gap:12px;align-items:center}.region-title{border:0;background:none;padding:0;font-size:18px;font-weight:720;text-align:left}.source-desc{font-size:12px;color:var(--muted);margin-top:3px}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:8px;margin:14px 0 3px}.metrics strong{display:block;font-size:22px;letter-spacing:-.5px;font-weight:700}.metrics span{font-size:11px;color:var(--muted)}.workbench{display:grid;grid-template-columns:minmax(0,1.8fr) minmax(320px,1fr);gap:18px;align-items:start}.explorer{padding:18px}.inspector{padding:18px;position:sticky;top:16px}.toolrow{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:12px 0}.toolrow>label{font-size:12px;color:var(--muted);display:flex;gap:6px;align-items:center}button,input,select{font:inherit}button{cursor:pointer;border:1px solid #c9d4dc;border-radius:6px;background:white;color:var(--ink);padding:6px 11px;font-size:12px;line-height:1.45}button:hover{background:#edf3f6;border-color:#8295a3}button:focus-visible,input:focus-visible,select:focus-visible{outline:3px solid #edbd57;outline-offset:2px}button.primary{background:var(--ink);color:white;border-color:var(--ink)}button.primary:hover{background:#305268}button[aria-pressed="true"]{background:#e8f0f4;border-color:#789da8;font-weight:700}input[type=number]{width:91px;border:1px solid #c9d4dc;border-radius:5px;padding:5px 7px;font-size:12px}input[type=search]{min-width:150px;width:100%;max-width:270px;border:1px solid #c9d4dc;border-radius:5px;padding:7px 9px;font-size:12px}select{font-size:12px;padding:6px;border:1px solid #c9d4dc;border-radius:5px;background:white}.spacer{flex:1}.position-label{font-size:13px;color:var(--muted);margin:8px 0 0}.svg-shell{overflow-x:auto;min-width:0}.svg-shell svg{display:block;width:100%;height:auto;min-width:440px}.detail-shell svg{min-width:590px}.legend{display:flex;gap:15px;align-items:center;flex-wrap:wrap;font-size:11px;color:var(--muted);margin:6px 0 10px}.swatch{display:inline-block;width:15px;height:9px;border-radius:2px;vertical-align:middle;margin-right:4px}.swatch.selected{background:#fff3d2;border:2px solid #a06c00}.swatch.repeated{background:repeating-linear-gradient(45deg,#75889a 0px,#75889a 2px,#fff 2px,#fff 4px)}.divider{height:1px;background:var(--line);margin:17px 0}.tablewrap{overflow:auto;max-height:460px;border:1px solid var(--line);border-radius:6px}table{width:100%;border-collapse:collapse;font-size:11px;white-space:nowrap}th{position:sticky;top:0;z-index:1;background:#eef3f6;text-align:left;font-weight:700;padding:9px 9px;border-bottom:1px solid var(--line)}td{padding:8px 9px;border-bottom:1px solid #edf1f3}tbody tr{cursor:pointer}tbody tr:hover{background:#f0f6f8}tbody tr.selected{background:#fff3d5}td.coords{max-width:190px;overflow:hidden;text-overflow:ellipsis}.tag{font-size:10px;border-radius:4px;background:#ecf1f4;padding:2px 5px;font-weight:650}.right{text-align:right}.mono{font-family:ui-monospace,SFMono-Regular,Consolas,"Liberation Mono",monospace}.probe-name{font-size:18px;line-height:1.3;margin:6px 0;overflow-wrap:anywhere}.facts{display:grid;grid-template-columns:1fr 1fr;gap:8px 14px;margin:12px 0}.facts div{font-size:12px;color:var(--muted)}.facts strong{display:block;color:var(--ink);font-size:15px}.sequence{padding:11px;border-radius:6px;background:#f1f5f7;font:12px/1.9 ui-monospace,SFMono-Regular,Consolas,monospace;word-break:break-all;letter-spacing:.3px}.docking{color:#73808c}.binding-arm{color:#102f44;font-weight:750}.sequence-key{font-size:10px;color:var(--muted);margin-top:4px}.copy-list{display:flex;flex-wrap:wrap;gap:6px;max-height:200px;overflow:auto}.copy-list button{font-size:11px;padding:4px 7px}.copy-list button.active{background:#fff1c9;border-color:#ad841f}.duplex{overflow-x:auto;background:#f6f8fa;border:1px solid #e4eaf0;border-radius:6px;padding:10px 11px;margin-top:9px}.duplex pre{font:11px/1.7 ui-monospace,SFMono-Regular,Consolas,monospace;margin:0;white-space:pre}.strand-note{font-size:12px;background:#edf3f6;border-radius:6px;padding:10px 11px;line-height:1.6}.small-link{font-size:11px}.status{font-size:12px;min-height:18px;margin:3px 0;color:var(--muted)}.status.error{color:#a73825}details{font-size:12px;color:var(--muted)}details summary{cursor:pointer;color:var(--ink);font-weight:650}details p{margin:8px 0}.gap-list{display:flex;flex-wrap:wrap;gap:5px;margin:10px 0}.footer{margin-top:22px;padding:18px 20px}.footer-grid{display:grid;grid-template-columns:1fr 1fr;gap:28px}.footer h3{margin-top:0}.footer p{font-size:12px;color:var(--muted);line-height:1.65}.empty{padding:25px;text-align:center;color:var(--muted)}.range-info{font-size:11px;color:var(--muted)}.hidden{display:none!important}.over-map-help{font-size:10px;color:var(--muted);margin-top:2px}.noscript{padding:18px;background:#fff0d5;color:#4b3215}
@media(max-width:1100px){main{padding:24px 18px}.workbench{grid-template-columns:minmax(0,1.5fr) minmax(300px,1fr)}.metrics strong{font-size:19px}.metrics span{font-size:10px}.headrow{display:block}.badges{justify-content:flex-start}}
@media(max-width:880px){.workbench{grid-template-columns:1fr}.inspector{position:static}.overview{grid-template-columns:1fr}.footer-grid{grid-template-columns:1fr;gap:12px}h1{font-size:27px}.metrics span{font-size:11px}.metrics strong{font-size:22px}}
@media print{body{background:white}main{padding:0}.workbench{display:block}.inspector{position:static;margin-top:20px}.toolrow,button,.tablewrap{display:none}.card{break-inside:avoid}.overview{grid-template-columns:1fr 1fr}.svg-shell svg{min-width:0}.footer{break-before:page}}
</style>
</head>
<body>
<noscript><div class="noscript">This viewer requires JavaScript in your browser. All data are embedded; no network connection or installation is needed.</div></noscript>
<main>
<header>
  <div class="headrow"><div><div class="eyebrow">RASER-FISH / SEQUENCE MAP</div><h1>Primary probes / 5-base minimum gaps</h1></div><div class="badges" id="badges"></div></div>
  <p class="lede">Explore the redesigned DP7 and DP3 pools. Every exact repeat copy is included in the spacing check. Click a site to inspect its binding arm, target strand and all repeated copies.</p>
  <div class="notice"><strong>What is mapped:</strong> full-length, exact matches of the DNA-binding arms only. Coordinates start at 1 in each supplied sequence. These are predicted sequence sites, not measurements of hybridization or fluorescence.</div>
</header>
<div class="strand-note" style="margin-bottom:18px;font-size:14px" id="strandBanner"></div>
<div id="warnings" class="notice hidden"></div>
<section class="overview" id="overview" aria-label="Whole-sequence binding maps"></section>
<div class="workbench">
<section class="card explorer">
  <div class="cardheading"><h2>Explore a sequence</h2><div id="regionTabs" class="toolrow" style="margin:0"></div></div>
  <div class="toolrow"><label>Start <input type="number" id="startInput" min="1" aria-label="Window start coordinate"></label><label>End <input type="number" id="endInput" min="1" aria-label="Window end coordinate"></label><button class="primary" id="goButton">Go</button><button id="fullButton">Full sequence</button><button id="zoomIn" title="Zoom in twofold">Zoom +</button><button id="zoomOut" title="Zoom out twofold">Zoom -</button><button id="panLeft" aria-label="Move window left">&larr;</button><button id="panRight" aria-label="Move window right">&rarr;</button></div>
  <div class="position-label" id="positionLabel"></div>
  <div class="note" id="spacingSummary" style="margin:4px 0 10px"></div>
  <div id="detailMap" class="svg-shell detail-shell" aria-label="Zoomable binding-site map"></div>
  <div class="legend"><span><i class="swatch" style="background:var(--dp7)"></i>DP7</span><span><i class="swatch" style="background:var(--dp3)"></i>DP3</span><span><i class="swatch repeated"></i>Multiple exact sites in this sequence</span><span><i class="swatch selected"></i>Selected probe</span></div>
  <div class="toolrow"><label><input type="checkbox" id="shadeGaps" checked> Shade unprobed spans</label><button id="showFiveGap">Inspect a 5-base gap</button><span class="spacer"></span><button id="exportSvg">Save this map (SVG)</button><button id="exportRegion">Export region sites (CSV)</button></div>
  <div id="status" class="status" role="status" aria-live="polite"></div>
  <details><summary id="gapSummary">Uncovered sequence</summary><p>Intervals with no full-arm exact site from any loaded probe. This does not imply the DNA is physically inaccessible.</p><div id="gapList" class="gap-list"></div></details>
  <div class="divider"></div>
  <div class="cardheading"><h2>Probe list</h2><span class="range-info" id="tableCount"></span></div>
  <div class="toolrow"><input type="search" id="searchInput" placeholder="Find probe name or arm sequence" aria-label="Search probe list"><select id="readoutFilter" aria-label="Filter probe list by readout"><option value="">All readouts</option></select><label><input type="checkbox" id="repeatsOnly"> Repeats only</label></div>
  <p class="range-info">One row per oligo. Filters apply to this list; maps continue to show every site.</p>
  <div class="tablewrap"><table><thead><tr><th>Probe</th><th>Readout</th><th class="right">Arm nt</th><th class="right">GC %</th><th class="right">Sites</th><th>Positions (inclusive)</th></tr></thead><tbody id="probeRows"></tbody></table></div>
</section>
<aside class="card inspector" id="inspector" aria-label="Selected probe details">
  <div class="eyebrow">PROBE INSPECTOR</div>
  <div id="inspectorContent"></div>
</aside>
</div>
<section class="card footer"><div class="footer-grid"><div><h3>Reading the strand labels</h3><p>The supplied 5&prime;-to-3&prime; sequence is the <strong>TOP (+)</strong> strand. The complementary 3&prime;-to-5&prime; row is the <strong>BOTTOM (-)</strong> strand. An arm identical to a window of that sequence binds the <strong>BOTTOM (-)</strong> DNA strand. An arm equal to the reverse complement of a window binds the supplied <strong>+</strong> strand. Both possibilities were searched independently.</p><p>Selected probe copies are highlighted in every overview. Probe bars are annotations below the DNA; their vertical offset is not a biological distance. Colored footprints on the DNA line indicate the strand bound. Overlapping inputs, if loaded, use extra lanes. Repeat status means multiple full-arm exact sites within a supplied sequence, not a genome-wide repeat annotation.</p></div><div><h3>Scope and export</h3><p>Each reference is mapped separately as a linear molecule. Their relative order, spacing, genomic integration and transcription boundaries are not inferred. No mismatches, short partial matches, circular-junction sites, RNA annotations, or genome-wide off-targets are included. Coverage is the union of matched intervals; overlapping bases are counted once.</p><p id="inputInfo"></p><div class="toolrow"><button id="exportAll">Export all binding sites (CSV)</button><button id="exportData">Save mapping data (JSON)</button></div></div></div></section>
</main>
<script id="dataset" type="application/json">__DATA_JSON__</script>
<script>
'use strict';
const DATA=JSON.parse(document.getElementById('dataset').textContent);
const refs=new Map(DATA.references.map(r=>[r.id,r]));
const palette={DP7:'#176e7b',DP3:'#7654ac',Other:'#687a8d'};
const state={region:DATA.references[0].id,start:1,end:DATA.references[0].length,selected:null,selectedHit:null,shadeGaps:true};
const $=id=>document.getElementById(id);
const fmt=n=>Number(n).toLocaleString('en-US');
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const comp=s=>s.replace(/[ACGTRYSWKMBDHVN]/g,c=>({A:'T',C:'G',G:'C',T:'A',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',N:'N'}[c]));
const color=p=>palette[p.readout]||palette.Other;
const regionHits=r=>r.hits.map(id=>DATA.hits[id]);
const probeHits=p=>p.hits.map(id=>DATA.hits[id]);
const countHere=(pi,rid)=>DATA.probes[pi].hits.filter(id=>DATA.hits[id].reference_id===rid).length;
const codepointSafeId=s=>Array.from(String(s)).map(c=>c.codePointAt(0).toString(16)).join('_');
const clamp=(n,min,max)=>Math.max(min,Math.min(max,n));
function status(text,error=false){$('status').textContent=text;$('status').className='status'+(error?' error':'');}
function setWindow(start,end){
  const r=refs.get(state.region),span=clamp(Math.round(end-start+1),1,r.length);
  start=clamp(Math.round(start),1,r.length-span+1);state.start=start;state.end=start+span-1;
}
function changeRegion(id){
  state.region=id;state.start=1;state.end=refs.get(id).length;
  const selected=state.selected===null?null:DATA.probes[state.selected];
  const first=selected&&probeHits(selected).find(h=>h.reference_id===id);
  if(first){state.selectedHit=first.id;}
  else{state.selected=null;state.selectedHit=null;}
  status('');renderAll();
}
function selectHit(id,zoom=true){
  const h=DATA.hits[id];state.selected=h.probe_index;state.selectedHit=id;
  if(state.region!==h.reference_id){state.region=h.reference_id;state.start=1;state.end=refs.get(h.reference_id).length;}
  if(zoom){const span=Math.max(160,h.end-h.start+1+90);setWindow((h.start+h.end-span+1)/2,(h.start+h.end+span-1)/2);}
  status('');renderAll();
}
function selectProbe(index){
  const p=DATA.probes[index],all=probeHits(p),inRegion=all.filter(h=>h.reference_id===state.region);
  const h=inRegion.find(h=>h.start<=state.end&&h.end>=state.start)||inRegion[0]||all[0];
  if(h)selectHit(h.id,true);else{state.selected=index;state.selectedHit=null;renderAll();}
}
function tickValues(start,end,maxTicks){
  const span=end-start+1,rough=span/maxTicks,magnitude=Math.pow(10,Math.floor(Math.log10(Math.max(1,rough))));
  const step=[1,2,5,10].map(n=>n*magnitude).find(n=>n>=rough)||10*magnitude;
  const vals=[start],edgeGap=span/(maxTicks+1);
  for(let t=Math.ceil((start+1)/step)*step;t<end;t+=step){
    if(t-start>=edgeGap&&end-t>=edgeGap)vals.push(Math.round(t));
  }
  if(end!==start)vals.push(end);return vals;
}
function mapMarkup(ref,start,end,mini=false){
  const width=mini?1000:1100,left=mini?185:245,right=30,plotW=width-left-right;
  const span=end-start+1,xBoundary=b=>left+(b-(start-1))/span*plotW;
  const xBase=b=>xBoundary(b-.5);
  const hits=regionHits(ref).filter(h=>h.start<=end&&h.end>=start);
  const groups=new Map();
  for(const h of hits){const p=DATA.probes[h.probe_index],key=p.readout+'|'+h.bound_strand;if(!groups.has(key))groups.set(key,[]);groups.get(key).push(h);}
  const prefix=(mini?'mini-':'detail-')+codepointSafeId(ref.id),showBases=!mini&&plotW/span>=9;
  const top=mini?43:(showBases?116:102),rowH=mini?17:23,barH=mini?11:15;
  let cursorY=top,groupsData=[];
  for(const [key,gh] of [...groups.entries()].sort()){
    gh.sort((a,b)=>a.start-b.start||a.end-b.end);
    const minLanes=1,laneEnds=Array(minLanes).fill(-Infinity),placed=[];
    for(const h of gh){let lane=0;for(let i=1;i<laneEnds.length;i++)if(laneEnds[i]<laneEnds[lane])lane=i;
      if(laneEnds[lane]>=h.start){lane=laneEnds.length;laneEnds.push(-Infinity);}
      laneEnds[lane]=h.end;placed.push({h,lane});}
    groupsData.push({key,placed,y:cursorY,lanes:laneEnds.length});cursorY+=rowH*laneEnds.length+13;
  }
  if(!groupsData.length)cursorY+=35;
  const axisY=cursorY+6,height=axisY+34;
  let svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${esc(ref.id)} exact binding sites, positions ${start} through ${end}" data-ref="${esc(ref.id)}" data-mini="${mini}" data-left="${left}" data-plotw="${plotW}" data-start="${start}" data-end="${end}" style="font-family:system-ui,-apple-system,sans-serif"><title>${esc(ref.id)}: ${start}-${end}; full-arm exact sites</title><defs><pattern id="${prefix}-repeat" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="2" height="6" fill="white" opacity="0.5"/></pattern></defs><rect width="${width}" height="${height}" fill="white"/>`;
  if(state.shadeGaps){for(const g of ref.gaps){if(g.start>end||g.end<start)continue;const a=xBoundary(Math.max(g.start,start)-1),b=xBoundary(Math.min(g.end,end));svg+=`<rect x="${a}" y="20" width="${Math.max(0,b-a)}" height="${axisY-17}" fill="#e7ebef" opacity="0.65"/>`;}}
  if(mini){
    svg+=`<text x="0" y="25" font-size="13" fill="#607180">TOP +, 5' to 3'</text><line x1="${left}" y1="21" x2="${width-right}" y2="21" stroke="#9baab5" stroke-width="1.2"/><path d="M ${width-right-6} 17 L ${width-right} 21 L ${width-right-6} 25" fill="none" stroke="#9baab5"/>`;
  }else if(showBases){
    svg+=`<text x="0" y="36" font-size="13" fill="#607180">TOP + / 5' to 3'</text><text x="0" y="61" font-size="13" fill="#607180">BOTTOM - / 3' to 5'</text>`;
    for(let b=start;b<=end;b++){const c=ref.sequence[b-1],x=xBase(b);const selectedStrands=new Set(hits.filter(h=>h.probe_index===state.selected&&h.start<=b&&h.end>=b).map(h=>h.bound_strand));for(const strand of selectedStrands)svg+=`<rect x="${xBoundary(b-1)}" y="${strand==='-'?45:20}" width="${plotW/span}" height="23" fill="#fff1c8"/>`;svg+=`<text x="${x}" y="36" text-anchor="middle" font-family="monospace" font-size="14" fill="#172b3b">${c}</text><text x="${x}" y="61" text-anchor="middle" font-family="monospace" font-size="14" fill="#607180">${comp(c)}</text>`;}
  }else{
    for(const [y,tag,direction] of [[28,"TOP + / 5' to 3'",1],[54,"BOTTOM - / 3' to 5'",-1]]){svg+=`<text x="0" y="${y+4}" font-size="13" fill="#607180">${tag}</text><line x1="${left}" y1="${y}" x2="${width-right}" y2="${y}" stroke="#9baab5" stroke-width="1.2"/>`;const ax=direction===1?width-right:left,sign=direction===1?-1:1;svg+=`<path d="M ${ax+sign*6} ${y-4} L ${ax} ${y} L ${ax+sign*6} ${y+4}" fill="none" stroke="#9baab5"/>`;}
  }
  if(!mini){
    // Footprints on the physically bound DNA strand, not on the identity row.
    if(!showBases){for(const h of hits){const p=DATA.probes[h.probe_index],a=xBoundary(Math.max(start,h.start)-1),b=xBoundary(Math.min(end,h.end)),y=h.bound_strand==='-'?54:28;svg+=`<line x1="${a}" y1="${y}" x2="${b}" y2="${y}" stroke="${color(p)}" stroke-width="4" opacity="0.85"/>`;}}
    for(const g of (ref.adjacent_site_gaps||[])){
      if(g.length<0||g.start<start||g.end>end)continue;
      const a=xBoundary(g.start-1),b=xBoundary(g.end),gw=b-a,y=showBases?86:77;
      if(gw<24)continue;
      svg+=`<g data-gap-start="${g.start}" data-gap-end="${g.end}" style="cursor:pointer"><title>${g.length} unoccupied bases: ${g.start}-${g.end}</title><path d="M ${a} ${y-4} L ${a} ${y} L ${b} ${y} L ${b} ${y-4}" fill="none" stroke="${g.pass?'#577b5f':'#b6342e'}"/><text x="${(a+b)/2}" y="${y+13}" text-anchor="middle" font-size="11" fill="${g.pass?'#426c4b':'#b6342e'}">${g.length} bp</text></g>`;
    }
  }
  for(const group of groupsData){const [readout,strand]=group.key.split('|');svg+=`<text x="0" y="${group.y+barH+2}" font-size="${mini?13:14}" font-weight="650" fill="${palette[readout]||palette.Other}">${esc(readout)} binds ${strand==='-'?'BOTTOM (-)':'TOP (+)'}</text>`;
    for(const {h,lane} of group.placed){const p=DATA.probes[h.probe_index],a=xBoundary(Math.max(start,h.start)-1),b=xBoundary(Math.min(end,h.end)),bw=Math.max(1,b-a),y=group.y+lane*rowH;
      const selected=h.probe_index===state.selected,repeated=countHere(h.probe_index,ref.id)>1;
      const title=`${p.name}\n${ref.id}:${fmt(h.start)}-${fmt(h.end)}\nBinds ${h.bound_strand==='-'?'BOTTOM (-)':'TOP (+)'} strand | ${p.arm_length} nt arm\n${countHere(h.probe_index,ref.id)} exact site(s) in this sequence`;
      svg+=`<g class="probe-site" data-hit="${h.id}" tabindex="0" role="button" aria-label="${esc(title.replace(/\n/g,'; '))}" style="cursor:pointer"><title>${esc(title)}</title>`;
      if(selected)svg+=`<rect x="${a-2}" y="${y-3}" width="${bw+4}" height="${barH+6}" rx="2" fill="#fff2cc" stroke="#ae7c15" stroke-width="${mini?1.7:2.3}"/>`;
      svg+=`<rect x="${a}" y="${y}" width="${bw}" height="${barH}" rx="1" fill="${color(p)}"/>`;
      if(repeated)svg+=`<rect x="${a}" y="${y}" width="${bw}" height="${barH}" rx="1" fill="url(#${prefix}-repeat)"/>`;
      if(!mini&&bw>70){const text=p.name.replace(/^R\d+_DP\d+_(?:OS|G\d+)_/, '');svg+=`<text x="${a+bw/2}" y="${y+barH-3}" fill="white" text-anchor="middle" font-size="11" font-weight="700" pointer-events="none">${esc(text)}</text>`;}
      svg+='</g>';
    }
  }
  if(!hits.length)svg+=`<text x="${left}" y="${top+21}" font-size="14" fill="#607180">No full-arm exact sites in this window.</text>`;
  svg+=`<line x1="${left}" y1="${axisY}" x2="${width-right}" y2="${axisY}" stroke="#97a5b1"/>`;
  for(const t of tickValues(start,end,mini?6:8)){const x=xBase(t);svg+=`<line x1="${x}" y1="${axisY}" x2="${x}" y2="${axisY+5}" stroke="#97a5b1"/><text x="${x}" y="${axisY+22}" font-size="${mini?13:14}" text-anchor="${t===start?'start':t===end?'end':'middle'}" fill="#607180">${fmt(t)}</text>`;}
  svg+=`<text x="0" y="${axisY+22}" font-size="${mini?12:13}" fill="#607180">Position (1-based)</text>`;
  if(mini&&state.region===ref.id&&(state.start>1||state.end<ref.length)){const a=xBoundary(state.start-1),b=xBoundary(state.end);svg+=`<rect x="${a}" y="7" width="${Math.max(2,b-a)}" height="${axisY+3}" fill="#b1ccd6" fill-opacity="0.13" stroke="#697f8d" stroke-width="1.5" stroke-dasharray="4 4" pointer-events="none"/>`;}
  return svg+'</svg>';
}
function renderOverview(){
  $('overview').innerHTML=DATA.references.map(r=>{const s=r.summary,active=r.id===state.region;return `<article class="card overview-card${active?' active':''}"><div class="cardheading"><div><button class="region-title" data-region="${esc(r.id)}">${esc(r.id)} <span class="note">/ ${fmt(r.length)} nt</span></button><div class="source-desc">${esc(r.description.replace(/^\S+\s*/,''))}</div></div><span class="badge">${s.repeat_probe_rows} repeat oligos</span></div><div class="metrics"><div><strong>${s.unique_probe_rows}</strong><span>primary oligos</span></div><div><strong>${s.exact_sites}</strong><span>exact binding sites</span></div><div><strong>${s.coverage_percent.toFixed(1)}%</strong><span>sequence covered</span></div><div><strong>${fmt(s.covered_bases)}</strong><span>bases covered</span></div></div><div class="svg-shell">${mapMarkup(r,1,r.length,true)}</div><div class="over-map-help">Click a probe to inspect it. Click empty map space to center a 500-nt window.</div></article>`;}).join('');
}
function renderControls(){
  const r=refs.get(state.region);$('regionTabs').innerHTML=DATA.references.map(x=>`<button data-region="${esc(x.id)}" aria-pressed="${x.id===r.id}">${esc(x.id)}</button>`).join('');
  $('startInput').max=r.length;$('endInput').max=r.length;$('startInput').value=state.start;$('endInput').value=state.end;
  const n=regionHits(r).filter(h=>h.start<=state.end&&h.end>=state.start).length;
  $('positionLabel').textContent=`${r.id}: ${fmt(state.start)}-${fmt(state.end)} | ${fmt(state.end-state.start+1)} nt window | ${n} exact sites intersect this window`;
  $('detailMap').innerHTML=mapMarkup(r,state.start,state.end,false);
  $('spacingSummary').textContent=`All mapped copies: ${r.summary.spacing_pairs_checked} neighbor pairs checked; minimum gap ${r.summary.minimum_gap_nt??'n/a'} bp; ${r.summary.spacing_violations} violations. Gap = next start - previous end - 1.`;
  $('gapSummary').textContent=`Uncovered sequence: ${fmt(r.length-r.summary.covered_bases)} bases; longest gap ${fmt(r.summary.longest_uncovered_interval)} nt`;
  const gaps=[...r.gaps].sort((a,b)=>b.length-a.length).slice(0,20);
  $('gapList').innerHTML=gaps.map(g=>`<button data-gap-start="${g.start}" data-gap-end="${g.end}">${fmt(g.start)}-${fmt(g.end)} (${g.length} nt)</button>`).join('')+(r.gaps.length>20?'<span class="note">20 longest shown. All intervals are in the bundled CSV.</span>':'');
}
function renderTable(){
  const query=$('searchInput').value.trim().toUpperCase(),readout=$('readoutFilter').value,repeats=$('repeatsOnly').checked;
  const listed=DATA.probes.filter(p=>{const n=countHere(p.index,state.region);return n>0&&(!readout||readout===p.readout)&&(!repeats||n>1)&&(!query||p.name.toUpperCase().includes(query)||p.arm.includes(query));});
  $('tableCount').textContent=`${listed.length} / ${refs.get(state.region).summary.unique_probe_rows} oligos`;
  $('probeRows').innerHTML=listed.map(p=>{const hs=probeHits(p).filter(h=>h.reference_id===state.region).sort((a,b)=>a.start-b.start||a.bound_strand.localeCompare(b.bound_strand));const coords=hs.map(h=>`${fmt(h.start)}-${fmt(h.end)} (${h.bound_strand})`).join('; ');return `<tr data-probe="${p.index}" tabindex="0" role="button" aria-label="Inspect ${esc(p.name)}" class="${state.selected===p.index?'selected':''}"><td class="mono">${esc(p.name)}</td><td><span class="tag" style="color:${color(p)}">${esc(p.readout)}</span></td><td class="right">${p.arm_length}</td><td class="right">${p.gc_percent.toFixed(1)}</td><td class="right"><strong>${hs.length}</strong></td><td class="coords mono" title="${esc(coords)}">${esc(coords)}</td></tr>`;}).join('')||'<tr><td colspan="6" class="empty">No probes match these filters.</td></tr>';
}
function renderInspector(){
  if(state.selected===null){$('inspectorContent').innerHTML='<h2 style="margin-top:10px">Select a probe</h2><p class="note">Click a bar on either map, or a probe row in the list. All exact copies of that oligo will be highlighted.</p>';return;}
  const p=DATA.probes[state.selected],all=probeHits(p).sort((a,b)=>a.reference_id.localeCompare(b.reference_id)||a.start-b.start||a.bound_strand.localeCompare(b.bound_strand)),h=state.selectedHit===null?null:DATA.hits[state.selectedHit];
  const nHere=all.filter(h=>h.reference_id===state.region).length,plus=all.filter(h=>h.bound_strand==='+').length,minus=all.filter(h=>h.bound_strand==='-').length;
  let content=`<h2 class="probe-name mono">${esc(p.name)}</h2><span class="tag" style="color:${color(p)}">${esc(p.readout)}</span> <span class="note">${esc(p.pool||'Loaded probe')}</span><div class="facts"><div>Binding-arm length<strong>${p.arm_length} nt</strong></div><div>Arm GC content<strong>${p.gc_percent.toFixed(1)}%</strong></div><div>Exact sites in ${esc(state.region)}<strong>${nHere}</strong></div><div>All supplied sequences<strong>${all.length} sites</strong></div></div><div class="strand-note"><strong>DNA strand bound: ${minus&&plus?'both + and -':minus?'BOTTOM / opposite (-)':plus?'TOP / supplied (+)':'no exact match'}</strong><br>${minus&&!plus?'The arm is identical to the supplied forward window, so it pairs with the BOTTOM, complementary DNA strand (3\' to 5\' left-to-right).':plus&&!minus?'The arm is the reverse complement of the supplied forward window, so it pairs with the TOP, supplied DNA strand.':'Both orientations were searched separately. See the strand label on each individual site.'}</div>`;
  content+=`<h3>DNA-binding arm <span class="note">5' to 3'</span></h3><div class="sequence">${esc(p.arm)}</div><div class="toolrow"><button id="copyArm">Copy arm</button><button id="showCopies">Show all copies in region</button></div>`;
  content+=`<h3>Every exact binding site</h3><div class="copy-list">${all.map(site=>`<button data-inspect-hit="${site.id}" class="${site.id===state.selectedHit?'active':''}">${esc(site.reference_id)}:${fmt(site.start)}-${fmt(site.end)} (${site.bound_strand})</button>`).join('')||'<span class="note">No matches in the supplied references.</span>'}</div>`;
  if(h){const r=refs.get(h.reference_id),window=r.sequence.slice(h.start-1,h.end),target=h.bound_strand==='-'?comp(window):window,armAligned=h.bound_strand==='-'?p.arm:p.arm.split('').reverse().join(''),leftEnd=h.bound_strand==='-'?"3'":"5'",rightEnd=h.bound_strand==='-'?"5'":"3'",pLeft=h.bound_strand==='-'?"5'":"3'",pRight=h.bound_strand==='-'?"3'":"5'";
    const line=(label,left,sequence,right)=>label.padEnd(18,' ')+left+' '+sequence+' '+right;
    const alignment=line('TOP (+) reference',"5'",window,"3'")+"\n\n"+line(h.bound_strand==='-'?'BOTTOM (-) DNA':'TOP (+) DNA',leftEnd,target,rightEnd)+"\n"+' '.repeat(21)+'|'.repeat(p.arm_length)+"\n"+line('Primary arm',pLeft,armAligned,pRight);
    content+=`<h3>Spacing at this binding site</h3><div class="facts"><div>Gap to previous site<strong>${h.gap_before_nt===null?'First site':h.gap_before_nt+' bp'}</strong></div><div>Gap to next site<strong>${h.gap_after_nt===null?'Last site':h.gap_after_nt+' bp'}</strong></div></div><p class="note">Measured along DNA; includes neighboring copies of other probes and of this same probe. Terminal flanks are not inter-probe gaps.</p>`;
    content+=`<h3>Selected-site duplex</h3><div class="note">${esc(h.reference_id)}:${fmt(h.start)}-${fmt(h.end)}; both strands shown left-to-right in reference coordinates.</div><div class="duplex"><pre>${esc(alignment)}</pre></div>`;
  }
  if(p.full){content+=`<h3>Full primary <span class="note">${p.full_length} nt; 5' to 3'</span></h3><div class="sequence"><span class="docking">${esc(p.prefix)}</span><span class="binding-arm">${esc(p.arm)}</span><span class="docking">${esc(p.suffix)}</span></div><div class="sequence-key">Bold = DNA-binding arm. Grey = non-targeting flanks, including docking handles and TT spacers in the bundled probes.</div><div class="toolrow"><button id="copyFull">Copy full primary</button></div>`;}
  content+='<p class="note" style="margin-top:14px">The duplex is a sequence-based illustration. It does not establish DNA accessibility, actual occupancy or absence of RNA binding.</p>';
  $('inspectorContent').innerHTML=content;
}
function renderAll(){renderOverview();renderControls();renderTable();renderInspector();}
function download(name,text,type){const url=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),2000);}
const safeFilename=s=>String(s).replace(/[^A-Za-z0-9._-]/g,'_');
const csvCell=x=>{let t=String(x??'');if(t.length>1&&/^[=+@]/.test(t))t="'"+t;return '"'+t.replace(/"/g,'""')+'"';};
function exportSites(onlyRegion){
  const header=['Probe_Name','Readout','Reference','Start_1based','End_1based_inclusive','DNA_strand_bound','Arm_match_to_forward_reference','Arm_length_nt','Arm_5to3','Gap_before_nt','Gap_after_nt'];
  const hs=DATA.hits.filter(h=>!onlyRegion||h.reference_id===state.region).sort((a,b)=>a.reference_id.localeCompare(b.reference_id)||a.start-b.start||a.end-b.end);
  const rows=hs.map(h=>{const p=DATA.probes[h.probe_index];return[p.name,p.readout,h.reference_id,h.start,h.end,h.bound_strand,h.arm_match_to_forward,p.arm_length,p.arm,h.gap_before_nt,h.gap_after_nt];});
  download(onlyRegion?safeFilename(state.region)+'_binding_sites.csv':'all_binding_sites.csv',[header,...rows].map(row=>row.map(csvCell).join(',')).join('\r\n')+'\r\n','text/csv;charset=utf-8');
  status(`Exported ${hs.length} exact sites.`);
}
async function copyText(text){
  try{if(navigator.clipboard&&window.isSecureContext){await navigator.clipboard.writeText(text);}else{const a=document.createElement('textarea');a.value=text;a.style.position='fixed';a.style.opacity='0';document.body.append(a);a.select();const ok=document.execCommand('copy');a.remove();if(!ok)throw new Error('Copy unavailable');}status('Sequence copied to clipboard.');}
  catch(e){status('Automatic copying was blocked. Select the displayed sequence and copy it manually.',true);}
}
$('goButton').addEventListener('click',()=>{const s=Number($('startInput').value),e=Number($('endInput').value),r=refs.get(state.region);if(!Number.isInteger(s)||!Number.isInteger(e)||s<1||e<s||e>r.length){status(`Enter integers satisfying 1 <= start <= end <= ${fmt(r.length)}.`,true);return;}setWindow(s,e);status('');renderAll();});
for(const id of ['startInput','endInput'])$(id).addEventListener('keydown',e=>{if(e.key==='Enter')$('goButton').click();});
$('fullButton').addEventListener('click',()=>{state.start=1;state.end=refs.get(state.region).length;status('');renderAll();});
function zoom(factor){const r=refs.get(state.region),oldSpan=state.end-state.start+1,span=clamp(Math.round(oldSpan*factor),Math.min(30,r.length),r.length),center=(state.start+state.end)/2;setWindow(center-(span-1)/2,center+(span-1)/2);status('');renderAll();}
$('zoomIn').addEventListener('click',()=>zoom(.5));$('zoomOut').addEventListener('click',()=>zoom(2));
for(const [id,direction] of [['panLeft',-1],['panRight',1]])$(id).addEventListener('click',()=>{const move=Math.round((state.end-state.start+1)*.75)*direction;setWindow(state.start+move,state.end+move);status('');renderAll();});
$('shadeGaps').addEventListener('change',()=>{state.shadeGaps=$('shadeGaps').checked;renderAll();});
$('searchInput').addEventListener('input',renderTable);$('readoutFilter').addEventListener('change',renderTable);$('repeatsOnly').addEventListener('change',renderTable);
$('showFiveGap').addEventListener('click',()=>{const r=refs.get(state.region),g=(r.adjacent_site_gaps||[]).find(g=>g.length===DATA.required_min_gap_nt);if(!g){status('No gap exactly equal to the required minimum in this region.');return;}const h=DATA.hits[g.left_hit_id];state.selected=h.probe_index;state.selectedHit=h.id;setWindow(g.start-37,g.end+37);renderAll();status(`Gap ${g.start}-${g.end}: ${g.length} unoccupied DNA bases.`);});
$('exportRegion').addEventListener('click',()=>exportSites(true));$('exportAll').addEventListener('click',()=>exportSites(false));
$('exportData').addEventListener('click',()=>download('binding_map.json',JSON.stringify(DATA,null,2),'application/json'));
$('exportSvg').addEventListener('click',()=>{let svg=$('detailMap').querySelector('svg').cloneNode(true);svg.removeAttribute('data-ref');const text=new XMLSerializer().serializeToString(svg);download(safeFilename(state.region)+`_${state.start}-${state.end}_probe_map.svg`,text,'image/svg+xml');status('Saved the current sequence map as SVG.');});
document.addEventListener('click',event=>{
  const t=event.target,site=t.closest('[data-hit]'),ins=t.closest('[data-inspect-hit]'),region=t.closest('[data-region]'),probe=t.closest('[data-probe]'),gap=t.closest('[data-gap-start]');
  if(site){selectHit(Number(site.dataset.hit));return;}if(ins){selectHit(Number(ins.dataset.inspectHit));return;}if(region){changeRegion(region.dataset.region);return;}if(probe){selectProbe(Number(probe.dataset.probe));return;}
  if(gap){setWindow(Number(gap.dataset.gapStart)-40,Number(gap.dataset.gapEnd)+40);status('');renderAll();return;}
  if(t.id==='copyArm'&&state.selected!==null){copyText(DATA.probes[state.selected].arm);return;}if(t.id==='copyFull'&&state.selected!==null){copyText(DATA.probes[state.selected].full);return;}
  if(t.id==='showCopies'&&state.selected!==null){const hs=probeHits(DATA.probes[state.selected]).filter(h=>h.reference_id===state.region);if(hs.length){setWindow(Math.min(...hs.map(h=>h.start))-35,Math.max(...hs.map(h=>h.end))+35);renderAll();}return;}
  const map=t.closest('svg[data-mini="true"]');if(map){const rect=map.getBoundingClientRect(),vw=map.viewBox.baseVal.width,x=(event.clientX-rect.left)/rect.width*vw,left=Number(map.dataset.left),pw=Number(map.dataset.plotw);if(x<left||x>left+pw)return;const rid=map.dataset.ref,r=refs.get(rid),center=1+(x-left)/pw*(r.length-1);if(rid!==state.region){state.region=rid;state.selected=null;state.selectedHit=null;}const span=Math.min(500,r.length);setWindow(center-(span-1)/2,center+(span-1)/2);renderAll();}
});
document.addEventListener('keydown',event=>{const t=event.target;if((event.key==='Enter'||event.key===' ')&&(t.matches('[data-hit]')||t.matches('[data-probe]'))){event.preventDefault();t.dispatchEvent(new MouseEvent('click',{bubbles:true}));}});
$('badges').innerHTML=`<span class="badge">${DATA.probes.length} primary oligos</span><span class="badge">${DATA.hits.length} exact sites</span><span class="badge good">${DATA.required_min_gap_nt}-base gaps: ${DATA.spacing_all_pass?'PASS':'CHECK'}</span><span class="badge">Offline</span>`;
const allBottom=DATA.hits.length>0&&DATA.hits.every(h=>h.bound_strand==='-');
$('strandBanner').innerHTML=allBottom?'<strong>These probes bind the BOTTOM (-) DNA strand.</strong> The TOP (+) row is the supplied sequence, identical in polarity to the primary binding arms. The spacing result shown above is recomputed from all mapped sites.':'<strong>Strands are labeled TOP (+) and BOTTOM (-).</strong> See each site for its bound strand and measured gaps.';
$('readoutFilter').innerHTML='<option value="">All readouts</option>'+[...new Set(DATA.probes.map(p=>p.readout))].sort().map(x=>`<option value="${esc(x)}">${esc(x)}</option>`).join('');
$('inputInfo').textContent=`Probe input: ${DATA.probe_source_file||'supplied table'}. The Python source and input files are included in the complete package. Re-run the script to load a different probe or reference set.`;
if(DATA.warnings.length){$('warnings').classList.remove('hidden');$('warnings').textContent=DATA.warnings.join(' ');}
if(DATA.hits.length){state.selected=DATA.hits[0].probe_index;state.selectedHit=DATA.hits[0].id;}
renderAll();
// Exposed read-only-by-convention for local verification and browser-console inspection.
window.probeViewer={data:DATA,state:state,selectHit:selectHit,selectProbe:selectProbe,changeRegion:changeRegion,render:renderAll};
</script>
</body>
</html>

"""

if __name__ == "__main__":
    sys.exit(main())
