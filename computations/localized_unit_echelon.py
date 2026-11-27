"""Row reduction using only units of a specified principal open set."""
from collections import Counter


def unit_echelon(rows,units,record=None):
    typ = next((type(v) for row in rows for v in row.values()),None)
    if typ is None:
        return {}

    unit_cache = {}
    def nonunit_part(f):
        while f:
            common = f.gcd(units)
            if common.is_constant():
                return f/f.leading_coefficient()
            f = f//common
        return f
    def is_unit(f):
        if f.is_constant():
            return bool(f)
        key = str(f/f.leading_coefficient())
        if key not in unit_cache:
            unit_cache[key] = nonunit_part(f).is_constant()
        return unit_cache[key]

    def require_unit(f):
        if not is_unit(f):
            assert record is not None,'Pivot outside the specified principal open set'
            record(nonunit_part(f))

    def primitive(row):
        if not row:
            return row
        vals = sorted(row.values(),key=len)
        common = vals[0]
        for f in vals[1:]:
            if common.is_constant():
                break
            common = common.gcd(f)
        if not common.is_constant():
            divisor = common//nonunit_part(common)
            if not divisor.is_constant():
                row = {j:v//divisor for j,v in row.items()}
        scalar = next(iter(row.values())).leading_coefficient()
        return {j:v/scalar for j,v in row.items()} if scalar!=1 else row

    def clear(row,pivot,j):
        if j not in row:
            return row
        a,b = row[j],pivot[j]
        g = a.gcd(b)
        a,b = a//g,b//g
        out = {k:b*v for k,v in row.items() if k!=j}
        for k,v in pivot.items():
            if k==j:
                continue
            z = out.get(k,0)-a*v
            if z:
                out[k] = z
            else:
                out.pop(k,None)
        return primitive(out)

    work = []
    for row in rows:
        if not row:
            continue
        den = next(iter(row.values())).d**0
        for x in row.values():
            den = (den//den.gcd(x.d))*x.d
        require_unit(den)
        work.append(primitive({j:x.n*(den//x.d) for j,x in row.items()}))

    pivots = {}
    while work:
        columns = Counter(j for row in work for j in row)
        choices = [((0 if f.is_constant() else 1 if is_unit(f) else 2,
                     (len(row)-1)*(columns[j]-1),len(f),len(row),j),i,j)
                   for i,row in enumerate(work) for j,f in row.items() if record is not None or is_unit(f)]
        assert choices,'No pivot invertible throughout the specified open set'
        _,i,j = min(choices)
        pivot = work.pop(i)
        require_unit(pivot[j])
        pivots[j] = pivot
        work = [out for row in work if (out:=clear(row,pivot,j))]
    for j in reversed(pivots):
        for k in list(pivots[j]):
            if k!=j and k in pivots:
                pivots[j] = clear(pivots[j],pivots[k],k)
    for j,row in pivots.items():
        require_unit(row[j])
    return {j:{k:typ(v,row[j]) for k,v in row.items()} for j,row in pivots.items()}
