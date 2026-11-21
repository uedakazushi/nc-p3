#!/usr/bin/env python3
"""Verify the first Ext step and actual two-return components with a pure first step."""
import argparse
import json
from math import comb
from pathlib import Path

from frobenius_tensor import reconstructed_windows
from two_window_ext import CutExt, PathRow, add_scaled, pbw_overlap_counts
from verify_indexed_anticommutative_ore import window
from verify_nonkoszul_boundary import tensor_at
from verify_quartic_obstruction import exact_rank


def profile(cut, row, offset, endpoint):
    dimensions, maps = cut.complex(row, offset, endpoint)
    for first, second in zip(maps, maps[1:]):
        for column in first:
            composition = {}
            for i, c in column.items():
                add_scaled(composition, second[i], c)
            assert not composition
    ranks = [0, *(exact_rank(matrix) for matrix in maps), 0]
    ext = [d-ranks[i]-ranks[i+1] for i, d in enumerate(dimensions)]
    assert ext[3] == row.dimension(offset+4+endpoint)
    return ext


def run(windows, regular):
    cut = CutExt(*windows[:2])
    results = []
    for start in range(4):
        row = PathRow(cut.windows, start, 8)
        first = [profile(cut, row, -start, i) for i in range(4)]
        pure = all(not any(e[:3]) for e in first)
        if regular:
            assert pure and row.dims == [comb(d+3, 3) for d in range(9)]
        second = profile(cut, row, 4-start, start)
        if regular:
            assert second == [0, 0, 0, 165]
        results.append({"start": start, "first_Ext_by_endpoint": first,
                        "first_pure": pure, "second_top_module_Ext_at_degree_8": second})
        print(start, pure, second, flush=True)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    signed = [window([0]*9, [1, 0, 0], i) for i in range(4)]
    boundary = reconstructed_windows(tensor_at(0))
    report = {"signed_Ore": run(signed, True), "Boolean": run(boundary, False)}
    assert any(not r["first_pure"] for r in report["Boolean"])
    assert any(r["first_pure"] and any(r["second_top_module_Ext_at_degree_8"][:3])
               for r in report["Boolean"])
    quantum = [{(b, a): 1, (a, b): -(2+a+2*b)} for b in range(4) for a in range(b)]
    counts = {name: pbw_overlap_counts(windows) for name, windows in
              (("quantum_rank_12", [quantum]*4), ("signed_Ore", signed), ("Boolean", boundary))}
    assert counts["quantum_rank_12"] == counts["signed_Ore"] == [4]*4
    assert max(counts["Boolean"]) > 4
    report["PBW_leading_triple_counts"] = counts
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
