"""Sparse doublet windows: complete support list after three necessary growth tests."""
from pathlib import Path
from subprocess import check_output, run
from tempfile import TemporaryDirectory

from verify_klein_mixed_positive import bit
from verify_klein_single_missing_supports import admissible
from verify_klein_three_lines import WORDS


def candidates(direction):
    name = 'klein_doublet_'+('parallel' if direction==1 else 'transverse')+'_masks.txt'
    return tuple(int(m,16) for m in Path(__file__).with_name(name).read_text().splitlines())


def main():
    source = Path(__file__).with_name('klein_two_doublet_supports.cpp')
    with TemporaryDirectory() as temp:
        binary = Path(temp)/'shapes'
        run(['c++','-std=c++17','-O3',str(source),'-o',str(binary)],check=True)
        for q,n in ((1,2779),(2,2794)):
            masks = set(map(int,check_output([str(binary),str(q)],text=True).splitlines()))
            assert len(masks)==n
            def rank(p,b):
                return 2 if b in (0,q if p%2 else 1) else 1
            surviving = [m for m in sorted(masks) if admissible([w for w in WORDS if m&bit(w)],rank)]
            assert tuple(surviving)==candidates(q)
            assert len(surviving)==(148 if q==1 else 51)
    print('PASS: both sparse doublet sectors, 2779 and 2794 orbit checks, 148 and 51 candidate supports')


if __name__=='__main__':
    main()
