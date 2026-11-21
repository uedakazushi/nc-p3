#!/usr/bin/env python3
"""Universal dimensions of one-return bimodule Tor, not matrix exactness."""
import argparse,json
from math import comb
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    hilbert=(1,4,6,4);records=[]
    for i in range(4):
        for j in range(4):
            degree=4+j-i;dimensions={}
            for h in range(degree):
                r=i-j+h
                if not 0<=r<=3:continue
                count=sum(comb(p+3,3)*comb(h-p+3,3) for p in range(4-i) if 0<=h-p<=j)
                if count:dimensions[h]=hilbert[r]*count
            euler=sum((-1)**(degree-1-h)*value for h,value in dimensions.items())
            assert euler==comb(degree+3,3)
            records.append({'start':i,'end':j,'internal_degree':degree,'dimensions':dimensions,'top_signed_euler':euler})
    assert records[6]['dimensions']=={1:8,2:144,3:480,4:400}
    assert records[5]['dimensions']=={0:1,1:32,2:156,3:160}
    assert records[3]['dimensions']=={3:120,4:1040,5:2400,6:1600}
    for start in range(4):
        for length in range(1,8):
            counts=[sum((start+n)%4==edge for n in range(length)) for edge in range(4)]
            assert 1 in counts
    assert all(sum(n%4==edge for n in range(8))==2 for edge in range(4))
    result={'result':'PASS','arithmetic':'exact integer dimensions and edge counts','all16_vertex_pairs':records,
            'one_return_cut_covers_all_internal_degrees_through7':True,'degree8_requires_multiple_return_parts':True,
            'limits':['No differential matrices or their ranks are calculated.','This is a degree<=7 Koszul truncation criterion, not a full Koszul or AS certificate.']}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
