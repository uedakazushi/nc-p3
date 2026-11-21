#!/usr/bin/env python3
"""Exact differentials at all cuts of a signed Ore control and a basepoint boundary."""
import argparse
import json
from math import comb
from pathlib import Path

from frobenius_tensor import reconstructed_windows
from general_cut_tor import CutFrames
from quintic_field_analysis import combine
from verify_indexed_anticommutative_ore import window
from verify_nonkoszul_boundary import tensor_at
from verify_quartic_obstruction import exact_rank


def verify(windows):
    frames = CutFrames(windows)
    components = {}
    for origin in range(4):
        for i in range(4):
            for j in range(4):
                degree = 4+j-i
                dims, maps = frames.complex(origin, i, j)
                for h in sorted(maps)[1:]:
                    assert all(not combine(maps[h-1], col, 0) for col in maps[h])
                ranks = {h: exact_rank(cols) for h, cols in maps.items()}
                homology = {h+1: dim-ranks.get(h, 0)-ranks.get(h+1, 0)
                            for h, dim in dims.items()}
                assert all(value >= 0 for value in homology.values())
                assert sum((-1)**(degree-n)*v for n, v in homology.items()) == comb(degree+3, 3)
                profile = {n: v for n, v in homology.items() if v}
                key = (degree, (origin+i) % 4)
                component = components.setdefault(key, {"homology": profile, "presentations": []})
                assert component["homology"] == profile, (key, origin, i, j)
                component["presentations"].append({"cut_origin": origin, "ranks": list(ranks.values())})
    assert len(components) == 28
    return [{"degree": degree, "start": start, **record}
            for (degree, start), record in sorted(components.items())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    signed_ore = verify([window([0]*9, [1, 0, 0], i) for i in range(4)])
    boundary = verify(reconstructed_windows(tensor_at(0)))
    for row in signed_ore:
        d = row["degree"]
        assert row["homology"] == {d: comb(d+3, 3)}
    defects = {5: [7, 0, 0, 0], 6: [21, 0, 0, 21], 7: [41, 0, 41, 63]}
    for row in boundary:
        d = row["degree"]
        defect = defects.get(d, [0]*4)[row["start"]]
        expected = {d: comb(d+3, 3)+defect}
        if defect:
            expected[d-1] = defect
        assert row["homology"] == expected
    result = {"arithmetic": "exact rational kernel frames, direct zero compositions and ranks",
              "tor_indices": "homological n in Tor_E(n,d); ranks in ascending nonzero-source h",
              "signed_Ore_H0_G100": signed_ore, "Boolean_t0": boundary}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    for row in boundary:
        if row["degree"] >= 5:
            print(row["degree"], row["start"], row["homology"])
    print("PASS two fixtures, all four cuts and 16 endpoint pairs per cut")


if __name__ == "__main__":
    main()
