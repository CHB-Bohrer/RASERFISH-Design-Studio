#!/usr/bin/env python3
"""Read-only OligoMiner-v1-compatible melting-temperature/composition audit.

This add-on preserves every input probe. It does not perform alignment, LDA,
Jellyfish, NUPACK, probe redesign, or experimental validation. Import
calculate_tm() and screen_arm() to use the same model during candidate selection.

OligoMiner source reference: beliveau-lab/OligoMiner, blockParse.py v1.7.
The temperature model, not the legacy Python-2 executable, is reproduced here
using explicitly pinned Biopython parameters. See HOW_TO_RUN.txt for scope.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import sys
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

try:
    import Bio
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError as exc:
    raise SystemExit("Install the dependency first: python3 -m pip install -r requirements.txt") from exc

VERSION = "1.0.0"
SOURCES = {
    "blockParse_v1_7": "https://github.com/beliveau-lab/OligoMiner/blob/master/blockParse.py",
    "probeTm_v1_7": "https://github.com/beliveau-lab/OligoMiner/blob/master/probeTm.py",
    "OligoMiner_README": "https://github.com/beliveau-lab/OligoMiner/blob/master/README.md",
    "Biopython_1_86": "https://biopython.org/docs/1.86/api/Bio.SeqUtils.MeltingTemp.html",
    "LoopTrace_methods": "https://doi.org/10.1038/s41467-025-61689-y",
}


@dataclass(frozen=True)
class MeltingProfile:
    name: str
    nn_table: str = "DNA_NN3"
    sodium_mM: float = 390.0
    higher_strand_nM: float = 25.0
    lower_strand_nM: float = 25.0
    formamide_percent: float = 50.0
    formamide_factor: float = 0.65
    saltcorr: int = 5
    round_before_formamide: bool = True


PROFILES = {
    "oligominer1_defaults": MeltingProfile(name="oligominer1_defaults"),
    "oligominer1_2nM": MeltingProfile(
        name="oligominer1_2nM", higher_strand_nM=2.0, lower_strand_nM=0.0),
    "previous_NN_audit": MeltingProfile(
        name="previous_NN_audit", nn_table="DNA_NN4", higher_strand_nM=2.0,
        lower_strand_nM=0.0, formamide_factor=0.61, round_before_formamide=False),
}


@dataclass(frozen=True)
class ScreeningPolicy:
    """Custom high-coverage audit policy; NOT all OligoMiner default settings."""
    minimum_length: int = 36
    maximum_length: int = 90
    minimum_gc_percent: float = 20.0
    maximum_gc_percent: float = 80.0
    minimum_tm_C: float = 42.0
    maximum_tm_C: float = 50.0
    prohibited_homopolymer_length: int = 5


def clean_dna(value: str, description: str = "sequence") -> str:
    seq = "".join(str(value).split()).upper()
    if not seq or set(seq) - set("ACGT"):
        raise ValueError(f"{description}: expected nonempty unmodified A/C/G/T sequence")
    return seq


def calculate_tm(arm: str, profile: MeltingProfile) -> dict[str, float]:
    """Perfect DNA-DNA duplex Tm; arm only, no handles or RNA energetics.

    Uses the explicit DNA_NN table and entropy salt correction in OligoMiner
    blockParse v1.7. That script rounds raw duplex Tm to 0.01 C before the
    formamide correction. It does not apply a self-complementarity correction.
    """
    arm = clean_dna(arm, "DNA-binding arm")
    if len(arm) < 2:
        raise ValueError("DNA-binding arm must have at least two bases")
    if profile.nn_table not in {"DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4"}:
        raise ValueError("Unsupported DNA thermodynamic table")
    vals = [profile.sodium_mM, profile.higher_strand_nM, profile.lower_strand_nM,
            profile.formamide_percent, profile.formamide_factor]
    if not all(math.isfinite(v) for v in vals):
        raise ValueError("All thermodynamic settings must be finite")
    c1, c2 = sorted([profile.higher_strand_nM, profile.lower_strand_nM], reverse=True)
    if profile.sodium_mM <= 0 or c2 < 0 or c1 - c2 / 2 <= 0:
        raise ValueError("Sodium must be positive; strand concentrations must give c1-c2/2 > 0")
    if not 0 <= profile.formamide_percent <= 100 or profile.formamide_factor < 0:
        raise ValueError("Formamide must be 0-100 percent with nonnegative correction factor")
    raw = mt.Tm_NN(arm, nn_table=getattr(mt, profile.nn_table),
                  Na=profile.sodium_mM, K=0, Tris=0, Mg=0, dNTPs=0,
                  dnac1=c1, dnac2=c2, selfcomp=False, saltcorr=profile.saltcorr)
    raw_screen = float(f"{raw:.2f}") if profile.round_before_formamide else raw
    correction = profile.formamide_percent * profile.formamide_factor
    return {"raw_duplex_C": raw, "adjusted_unrounded_C": raw - correction,
            "adjusted_screen_C": raw_screen - correction}


def screen_arm(arm: str, profile: MeltingProfile,
               policy: ScreeningPolicy = ScreeningPolicy()) -> dict[str, Any]:
    """Reusable candidate-stage gate. Does NOT assess specificity or folding.

    The standalone audit reports this outcome; it never deletes a sequence.
    Tm bounds are STRICT to match blockParse.py's tmCheck method.
    """
    arm = clean_dna(arm)
    tm_values = calculate_tm(arm, profile)
    gc = 100.0 * (arm.count("G") + arm.count("C")) / len(arm)
    max_run = max(len(m.group(0)) for m in re.finditer(r"A+|C+|G+|T+", arm))
    tests = {
        "length_pass": policy.minimum_length <= len(arm) <= policy.maximum_length,
        "gc_pass": policy.minimum_gc_percent <= gc <= policy.maximum_gc_percent,
        "tm_pass": policy.minimum_tm_C < tm_values["adjusted_screen_C"] < policy.maximum_tm_C,
        "homopolymer_pass": max_run < policy.prohibited_homopolymer_length,
    }
    return {**tm_values, "length": len(arm), "GC_percent": gc,
            "longest_homopolymer": max_run, **tests,
            "sequence_screen_pass": all(tests.values())}


def read_probes(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise ValueError("CSV must have unique, nonempty column names")
        rows = list(reader)
    if not rows:
        raise ValueError("Probe CSV is empty")
    if "DNA_Binding_Arm_5to3" not in fields:
        raise ValueError("Use a viewer/design audit CSV with DNA_Binding_Arm_5to3; do not use the IDT order sheet")
    if not ({"Probe_Name", "Name"} & set(fields)):
        raise ValueError("CSV must contain Probe_Name or Name")
    names = set()
    for n, row in enumerate(rows, 2):
        if None in row or any(v is None for v in row.values()):
            raise ValueError(f"CSV row {n}: wrong number of fields")
        name = row.get("Probe_Name") or row.get("Name")
        if not name or name in names or any(ch.isspace() for ch in name):
            raise ValueError(f"CSV row {n}: probe names must be unique, nonempty, and contain no whitespace")
        names.add(name)
        arm = clean_dna(row["DNA_Binding_Arm_5to3"], f"{name} arm")
        full = row.get("Full_Primary_5to3", "")
        if full and arm not in clean_dna(full, f"{name} full primary"):
            raise ValueError(f"{name}: arm is not contained in the supplied full primary")
    return rows, fields


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run_audit(input_path: Path, out: Path, active: MeltingProfile,
              policy: ScreeningPolicy) -> dict[str, Any]:
    rows, fields = read_probes(input_path)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory is not empty. Use a new directory to preserve previous results.")
    input_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
    profiles = dict(PROFILES)
    if active.name not in profiles or active != profiles[active.name]:
        profiles["active_custom"] = active
    computed = []
    summaries = []
    for row in rows:
        arm = clean_dna(row["DNA_Binding_Arm_5to3"])
        active_result = screen_arm(arm, active, policy)
        result = dict(row)  # preserve all original fields and sequence strings
        result.update({
            "OMQC_Active_profile": active.name,
            "OMQC_Arm_GC_percent": round(active_result["GC_percent"], 6),
            "OMQC_Longest_arm_homopolymer": active_result["longest_homopolymer"],
            "OMQC_36to42nt_length_only": 36 <= len(arm) <= 42,
            "OMQC_Active_adjusted_Tm_C": round(active_result["adjusted_screen_C"], 6),
            "OMQC_Active_length_pass": active_result["length_pass"],
            "OMQC_Active_GC_pass": active_result["gc_pass"],
            "OMQC_Active_Tm_pass": active_result["tm_pass"],
            "OMQC_Active_homopolymer_pass": active_result["homopolymer_pass"],
            "OMQC_Active_sequence_screen_pass": active_result["sequence_screen_pass"],
            "OMQC_Legacy_length_GC_Tm_C": round(64.9 + 41 *
                (arm.count("G") + arm.count("C") - 16.4) / len(arm) - 30.5, 6),
            "OMQC_Genome_alignment": "NOT_RUN",
            "OMQC_LDA_specificity": "NOT_RUN",
            "OMQC_Genome_kmer_screen": "NOT_RUN",
            "OMQC_Folding_and_dimers": "NOT_RUN",
            "OMQC_Probe_action": "RETAINED_UNCHANGED_AUDIT_ONLY",
        })
        for label, profile in profiles.items():
            values = calculate_tm(arm, profile)
            result[f"OMQC_{label}_Tm_C"] = round(values["adjusted_screen_C"], 6)
            result[f"OMQC_{label}_42ltTm_lt50"] = 42 < values["adjusted_screen_C"] < 50
        if len(result) != len(fields) + len(set(result) - set(fields)):
            raise ValueError("Unexpected duplicate input/output column")
        if any(key.startswith("OMQC_") for key in fields):
            raise ValueError("Input already has OMQC columns. Use the original design CSV, not an earlier audit.")
        computed.append(result)
    regions = sorted({r.get("Region", "ALL") or "ALL" for r in rows})
    for region in regions:
        subset = [r for r in rows if (r.get("Region", "ALL") or "ALL") == region]
        for label, profile in profiles.items():
            values = [calculate_tm(r["DNA_Binding_Arm_5to3"], profile)["adjusted_screen_C"] for r in subset]
            summaries.append({"Region": region, "Profile": label, "Oligos": len(subset),
                "Tm_min_C": round(min(values), 6), "Tm_median_C": round(statistics.median(values), 6),
                "Tm_max_C": round(max(values), 6),
                "Number_42ltTm_lt50": sum(42 < t < 50 for t in values),
                "Number_ge50C": sum(t >= 50 for t in values),
                "Number_le42C": sum(t <= 42 for t in values)})
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "OligoMiner_temperature_audit.csv", computed, list(computed[0]))
    write_csv(out / "temperature_summary.csv", summaries, list(summaries[0]))
    with (out / "binding_arms.fastq").open("w") as fq, (out / "binding_arms_for_probeTm.tsv").open("w") as tsv:
        for row in rows:
            name = row.get("Probe_Name") or row["Name"]
            arm = clean_dna(row["DNA_Binding_Arm_5to3"])
            fq.write(f"@{name}\n{arm}\n+\n{'I' * len(arm)}\n")
            tsv.write(f"{name}\t{arm}\n")
    manifest = {
        "tool": "RASERFISH OligoMiner thermodynamic QC add-on", "version": VERSION,
        "mode": "AUDIT_ONLY_NO_PROBES_REMOVED_OR_CHANGED", "input": str(input_path),
        "input_sha256": input_hash, "input_rows": len(rows), "output_rows": len(computed),
        "active_profile": asdict(active), "active_sequence_policy": asdict(policy),
        "comparison_profiles": {k: asdict(v) for k, v in profiles.items()},
        "software": {"python": sys.version.split()[0], "biopython": Bio.__version__},
        "sources": SOURCES,
        "screening_status": {"on_target_DNA_duplex_temperature": "RUN",
            "GC_length_homopolymer_flags": "RUN",
            "genome_alignment_LDA_kmers": "NOT_RUN", "NUPACK_folding": "NOT_RUN",
            "spacing": "NOT_RECALCULATED_BY_THIS_ADDON", "transcriptome_screen": "NOT_RUN"},
        "scope_notes": [
            "OligoMiner-v1-compatible temperature calculation, not execution of its entire pipeline.",
            "Default active policy is custom: lengths 36-90 and Tm strictly 42-50; GC 20-80 and no 5-base homopolymers.",
            "The 25/25-nM default is a software comparison, not a claim about fixed-cell target concentration.",
            "The 2/0-nM profile is a custom trace-target assumption, not an established in-situ concentration measurement.",
            "Exact published LoopTrace command line has not been recovered or reproduced.",
            "50% formamide is modeled empirically; predicted temperatures are not measured in-situ values.",
            "FASTQ preserves original names; unmodified OligoMiner outputClean expects coordinate-style IDs, so use a mapped ID adapter.",
            "No passed flag here is a genome-specificity, strand, spacing, or experimental validation certificate.",
        ],
        "summary": summaries,
    }
    (out / "audit_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == input_hash
    assert all(all(new[k] == old[k] for k in fields) for new, old in zip(computed, rows))
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--probes", type=Path, required=True, help="CSV with explicit DNA_Binding_Arm_5to3 and Probe_Name/Name")
    ap.add_argument("--out", type=Path, required=True, help="New or empty output directory")
    ap.add_argument("--profile", choices=list(PROFILES), default="oligominer1_2nM")
    ap.add_argument("--sodium-mM", type=float)
    ap.add_argument("--probe-nM", type=float)
    ap.add_argument("--target-nM", type=float)
    ap.add_argument("--formamide-percent", type=float)
    ap.add_argument("--min-tm", type=float, default=42)
    ap.add_argument("--max-tm", type=float, default=50)
    ap.add_argument("--min-gc", type=float, default=20)
    ap.add_argument("--max-gc", type=float, default=80)
    ap.add_argument("--min-length", type=int, default=36)
    ap.add_argument("--max-length", type=int, default=90)
    ap.add_argument("--reject-homopolymer-at", type=int, default=5)
    args = ap.parse_args()
    policy = ScreeningPolicy(args.min_length, args.max_length, args.min_gc, args.max_gc,
                             args.min_tm, args.max_tm, args.reject_homopolymer_at)
    if not (2 <= args.min_length <= args.max_length and 0 <= args.min_gc <= args.max_gc <= 100
            and args.min_tm < args.max_tm and args.reject_homopolymer_at >= 2):
        ap.error("Invalid length, GC, Tm, or homopolymer bounds")
    updates = {key: val for key, val in {
        "sodium_mM": args.sodium_mM, "higher_strand_nM": args.probe_nM,
        "lower_strand_nM": args.target_nM, "formamide_percent": args.formamide_percent,
    }.items() if val is not None}
    active = PROFILES[args.profile]
    if updates:
        active = replace(active, name=args.profile + "_custom", **updates)
    try:
        result = run_audit(args.probes, args.out, active, policy)
    except (OSError, ValueError) as exc:
        ap.error(str(exc))
    print(f"Audited {result['input_rows']} probes; all retained unchanged.")
    for row in result["summary"]:
        print(f"{row['Region']} / {row['Profile']}: "
              f"{row['Tm_min_C']:.2f}-{row['Tm_max_C']:.2f} C; "
              f"{row['Number_42ltTm_lt50']}/{row['Oligos']} inside 42<Tm<50 C")
    print("Not run: genome/transcriptome alignment, LDA, genome k-mers, or folding.")
    print(f"Audit files: {args.out.resolve()}")


if __name__ == "__main__":
    main()
