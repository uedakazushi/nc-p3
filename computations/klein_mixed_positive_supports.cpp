// Complete three even rectangles into two rank-two and two rank-one odd blocks.
#include "klein_support_shape.hpp"

struct Choice { int matrix, missing; };
struct Rectangle { int matrix; Four key{}; };

int main() {
    std::array<Shape,65536> shapes;
    for (int m=0; m<65536; ++m) shapes[m]=classify(m);
    std::array<std::array<std::vector<Choice>,4096>,4> options;
    for (int e=0; e<4; ++e) for (int m=1; m<65536; ++m) {
        if (!shapes[m].rank[e<2 ? 2 : 1]) continue;
        int key=0, missing=0;
        for (int j=0; j<4; ++j) for (int l=0; l<4; ++l) {
            if (!(m & (1<<(4*j+l)))) continue;
            int b=j^l^e;
            if (b) key |= 1<<(4*(b-1)+j);
            else missing |= 1<<(4*j+l);
        }
        options[e][key].push_back({m,missing});
    }
    std::array<std::vector<Rectangle>,4> rects;
    for (int b=1; b<4; ++b) for (int m=1; m<65536; ++m) {
        if (!shapes[m].rank[1]) continue;
        Rectangle r{m};
        for (int i=0; i<4; ++i) for (int k=0; k<4; ++k)
            if (m & (1<<(4*i+k))) r.key[i^k^b] |= 1<<(4*(b-1)+(i^b));
        rects[b].push_back(r);
    }
    std::uint64_t total=0, kept=0, dense=0;
    for (const auto& a:rects[1]) for (const auto& b:rects[2]) for (const auto& c:rects[3]) {
        Four key{};
        for (int e=0; e<4; ++e) key[e]=a.key[e]|b.key[e]|c.key[e];
        Four even{0,a.matrix,b.matrix,c.matrix};
        for (auto& o0:options[0][key[0]]) for (auto& o1:options[1][key[1]])
        for (auto& o2:options[2][key[2]]) for (auto& o3:options[3][key[3]]) {
            even[0]=o0.missing|o1.missing|o2.missing|o3.missing;
            if (!shapes[even[0]].rank[3]) continue;
            ++total;
            Four odd{o0.matrix,o1.matrix,o2.matrix,o3.matrix};
            std::array<int,16> edges{};
            for (int p=0; p<4; ++p) for (int h=0; h<4; ++h) {
                auto& sh=shapes[p%2 ? odd[h] : even[h]];
                int free=p<2 ? sh.free_rows : sh.free_columns;
                for (int v=0; v<4; ++v) if (free & (1<<v))
                    edges[4*p+v] |= 1<<(4*((p+1)%4)+(v^h));
            }
            if (has_branching_cycle(edges)) continue;
            ++kept;
            if ((a.matrix&b.matrix&c.matrix&o2.matrix&o3.matrix)==65535) { ++dense; continue; }
            std::uint64_t mask=0;
            for (int t=0; t<4; ++t) mask |= std::uint64_t(even[t])<<(16*t);
            std::cout<<mask<<'\n';
        }
    }
    if (total!=109034088 || kept!=15000 || dense!=256)
        throw std::logic_error("Unexpected mixed positive support counts");
}
