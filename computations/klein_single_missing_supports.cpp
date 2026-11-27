// Enumerate simultaneous rank supports and print masks with no exponential path subgraph.
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

struct Slice {
    int matrix = 0;
    std::array<int, 3> even_bits{};
};

std::array<Shape, 65536> shapes;
std::array<std::array<Slice, 4096>, 4> odd_slices;
std::array<std::array<std::vector<int>, 4096>, 3> completions;

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

int even_matrix(int b, int bits) {
    int result = 0;
    for (int e=0; e<4; ++e) for (int i=0; i<4; ++i)
        if (bits & (1 << (4*e+i))) result |= 1 << (4*i+(i^b^e));
    return result;
}

bool exponential(const std::array<int,3>& even, const Four& odd) {
    std::array<int,16> edges{};
    for (int i=0; i<4; ++i) {
        edges[i] = 1 << (4+i);
        edges[8+i] = 1 << (12+i);
    }
    for (int p=0; p<4; ++p) for (int c=(p%2 ? 0 : 1); c<4; ++c) {
        const Shape& sh = shapes[p%2 ? odd[c] : even[c-1]];
        int free = p<2 ? sh.free_rows : sh.free_columns;
        for (int i=0; i<4; ++i)
            if (free & (1<<i)) edges[p*4+i] |= 1 << (((p+1)%4)*4+(i^c));
    }
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

void enumerate(const Four& ranks) {
    std::array<std::vector<int>,4> choices;
    for (int e=0; e<4; ++e) for (int bits=0; bits<4096; ++bits)
        if (shapes[odd_slices[e][bits].matrix].rank[ranks[e]]) choices[e].push_back(bits);
    std::uint64_t count = 0, surviving = 0;
    for (int j:choices[1]) for (int k:choices[2]) for (int l:choices[3]) {
        std::array<int,3> partial{};
        for (int b=0; b<3; ++b)
            partial[b] = odd_slices[1][j].even_bits[b] + 16*odd_slices[2][k].even_bits[b] +
                         256*odd_slices[3][l].even_bits[b];
        for (int a:completions[0][partial[0]]) for (int b:completions[1][partial[1]])
            for (int c:completions[2][partial[2]]) {
                std::array<int,3> last{a,b,c}, even{};
                int m = 0;
                for (int t=0; t<3; ++t) {
                    for (int i=0; i<4; ++i)
                        if (last[t] & (1<<i)) m |= 1 << (4*t+(i^(t+1)));
                    even[t] = even_matrix(t+1,16*partial[t]+last[t]);
                }
                Four odd{odd_slices[0][m].matrix,odd_slices[1][j].matrix,
                         odd_slices[2][k].matrix,odd_slices[3][l].matrix};
                if (!shapes[odd[0]].rank[ranks[0]]) continue;
                ++count;
                if (exponential(even,odd)) continue;
                ++surviving;
                std::uint64_t code = m + (std::uint64_t(j)<<12) + (std::uint64_t(k)<<24) + (std::uint64_t(l)<<36);
                std::cout << code << '\n';
            }
    }
    if (ranks[0]==3 ? count!=15936576 || surviving!=0 : count!=115151588 || surviving!=2664)
        throw std::logic_error("Unexpected support count");
}

int main() {
    for (int mask=0; mask<65536; ++mask) shapes[mask] = classify(mask);
    for (int e=0; e<4; ++e) for (int bits=0; bits<4096; ++bits) {
        Slice& part = odd_slices[e][bits];
        for (int b=1; b<4; ++b) for (int j=0; j<4; ++j) {
            if (!(bits & (1 << (4*(b-1)+j)))) continue;
            part.matrix |= 1 << (4*j+(j^b^e));
            part.even_bits[b-1] |= 1 << (j^b);
        }
    }
    for (int b=1; b<4; ++b) for (int u=0; u<4096; ++u) for (int m=0; m<16; ++m)
        if (shapes[even_matrix(b,u*16+m)].rank[2]) completions[b-1][u].push_back(m);
    enumerate({3,1,1,1});
    enumerate({2,2,1,1});
}
