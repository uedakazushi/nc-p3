#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using Four = std::array<int, 4>;

struct Shape {
    std::array<bool, 4> rank{};
    int free_rows = 0;
    int free_columns = 0;
};

int root(int i, const std::array<int, 8>& parent) {
    while (parent[i] != i) i = parent[i];
    return i;
}

Shape classify(int mask) {
    Shape shape;
    Four rows{}, columns{};
    int active_rows = 0, active_columns = 0;
    for (int i=0; i<4; ++i) for (int j=0; j<4; ++j) {
        if (!(mask & (1 << (4*i+j)))) continue;
        rows[i] |= 1 << j;
        columns[j] |= 1 << i;
        active_rows |= 1 << i;
        active_columns |= 1 << j;
    }
    shape.free_rows = 15 ^ active_rows;
    shape.free_columns = 15 ^ active_columns;
    shape.rank[1] = mask && __builtin_popcount(mask) ==
        __builtin_popcount(active_rows)*__builtin_popcount(active_columns);

    std::array<int, 8> parent{0,1,2,3,4,5,6,7};
    for (int i=0; i<4; ++i) for (int j=0; j<4; ++j) {
        if (rows[i] && columns[j] && !(rows[i] & (1<<j)))
            parent[root(i,parent)] = root(4+j,parent);
    }
    bool consistent = true;
    int row_lines = 0, column_lines = 0;
    for (int i=0; i<4; ++i) {
        if (rows[i]) row_lines |= 1 << root(i,parent);
        if (columns[i]) column_lines |= 1 << root(4+i,parent);
        for (int j=0; j<4; ++j)
            if ((rows[i] & (1<<j)) && root(i,parent)==root(4+j,parent)) consistent = false;
    }
    shape.rank[2] = consistent && __builtin_popcount(row_lines)>=2 && __builtin_popcount(column_lines)>=2;
    Four permutation{0,1,2,3};
    int largest = 0, perfect = 0;
    do {
        int matches = 0;
        for (int i=0; i<4; ++i) matches += (rows[i] >> permutation[i]) & 1;
        largest = std::max(largest,matches);
        perfect += matches==4;
    } while (std::next_permutation(permutation.begin(),permutation.end()));
    shape.rank[3] = largest>=3 && perfect!=1;
    return shape;
}

bool has_branching_cycle(const std::array<int,16>& edges) {
    auto reach = edges;
    for (int k=0; k<16; ++k) for (int i=0; i<16; ++i)
        if (reach[i] & (1<<k)) reach[i] |= reach[k];
    for (int i=0; i<16; ++i) {
        int returning = 0;
        for (int targets=edges[i]; targets; targets &= targets-1)
            returning += (reach[__builtin_ctz(targets)] >> i) & 1;
        if (returning>1) return true;
    }
    return false;
}
