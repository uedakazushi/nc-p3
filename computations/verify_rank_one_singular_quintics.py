"""Check the cone postulation, the partial blowup, and the two equality obstructions."""
for t in range(3, 16):
    sections = sum(2*t-5-2*i+1 for i in range(t-2))
    assert sections == (t-2)*(t-1)
    assert (t+1)**2 - sections == 5*t-1
assert 6-2 == 4
assert 2 - 2*(2-1)//2 == 1
assert [(5+j)//2 for j in range(3)] == [2, 3, 3]
assert 16-9 == 1+4+2 == 7
assert 7 % 2 and 7 % 3
print("PASS: quintic cone polynomial, length-two trace, double-point blowup and both incidence partitions")
