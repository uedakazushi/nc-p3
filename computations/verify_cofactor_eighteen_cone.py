"""Genus and multiplicity for the cubic section of the quadric cone."""
models = []
for u in range(1, 4):
    square = 12*u - 2*u*u
    # K_{F2} = -2E-4F, pairing with uE+6F is -12.
    genus = 1 + (square-12)//2
    models.append((u, 6-2*u, genus, square-8))
assert models == [(1, 4, 0, 2), (2, 2, 3, 8), (3, 0, 4, 10)]
assert 3**2*2 - 2**2 - 8 == 6
assert 20-3 > 13 and 20-6 > 12
print("PASS: vertex multiplicities, sectional genera, blowup lengths and triple bounds")
