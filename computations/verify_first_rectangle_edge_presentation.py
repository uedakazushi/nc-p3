"""Regenerate polynomial certificates for freeness on the rank-ten support-85 surfaces."""
from pathlib import Path
import json
import sys

from flint_rational_function import rational_function_field
from localized_unit_echelon import unit_echelon
from quartic_center_presentation import certify_presentation
from verify_first_rectangle_edge_center import edge_rules
from verify_parallel_invertible_rectangles import s
from verify_rectangular_surface_cut_charts import relation_cut_data

CERTIFICATE = Path(__file__).resolve().parents[1]/'research'/'first_rectangle_edge_presentations.json'


def factors(start,alternate):
    v,z = s.symbols('v z')
    F = rational_function_field(('v','z'))
    units = F.from_sympy(v*z*(v-1)*(v+1)*(v+z)*(z-1)*(v*z+1)).n
    shifts = (2,0,2,0) if alternate else (0,2,2,0)
    constants = tuple(F.one*int(i==int(alternate)) for i in range(5))
    rel = edge_rules(v,z,F.from_sympy,shifts)
    _,_,_,sections = relation_cut_data(rel,F,units,full_center=True)
    excluded = set()
    def record(f):
        excluded.update(str(g) for g,_ in f.factor()[1])
    def pivot(rows):
        indices = []
        P = unit_echelon(rows,units,record,indices)
        return P,indices
    def span(rows):
        return unit_echelon(rows,units,record)
    certify_presentation(rel,sections,start,constants,F.one,pivot,span,(1,4,10,20,30,36,34,20,1,0))
    return sorted(excluded)


if __name__=='__main__':
    expected = json.loads(CERTIFICATE.read_text())
    for key in (sys.argv[1:] or expected):
        assert factors(int(key[-1]),key[0]=='B')==expected[key]
        print('PASS: all independent polynomial-module relations on chart',key,flush=True)
