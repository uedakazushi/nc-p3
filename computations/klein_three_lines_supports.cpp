// Simultaneous support tests for two adjacent (3,1,1,1) cuts.
#include "klein_support_shape.hpp"

struct Rectangle {
    int mask;
    Four crossing{};   // Row indices on each crossing block of the other cut.
    int thick = 0;
};

int main() {
    std::array<Shape,65536> shapes;
    for (int m=0; m<65536; ++m) shapes[m]=classify(m);
    std::array<std::vector<Rectangle>,4> rectangles;
    std::array<std::array<std::vector<int>,4096>,4> completions;
    for (int color=1; color<4; ++color) for (int m=1; m<65536; ++m) {
        if (!shapes[m].rank[1]) continue;
        Rectangle rect{m};
        int key=0;
        for (int i=0; i<4; ++i) for (int k=0; k<4; ++k) {
            if (!(m & (1<<(4*i+k)))) continue;
            int cross=i^k^color;
            if (cross) {
                key |= 1<<(4*(cross-1)+i);
                rect.crossing[cross] |= 1<<(i^cross);
            } else rect.thick |= 1<<(4*i+k);
        }
        rectangles[color].push_back(rect);
        completions[color][key].push_back(m);
    }
    std::uint64_t candidates=0, finite_paths=0, dense=0;
    for (const auto& a:rectangles[1]) for (const auto& b:rectangles[2])
    for (const auto& c:rectangles[3]) {
        Four key{};
        for (int t=1; t<4; ++t)
            key[t]=a.crossing[t]+16*b.crossing[t]+256*c.crossing[t];
        for (int i:completions[1][key[1]]) for (int j:completions[2][key[2]])
        for (int k:completions[3][key[3]]) {
            Four even{0,i,j,k}, odd{0,a.mask,b.mask,c.mask};
            int even_zero=a.thick|b.thick|c.thick, odd_zero=0;
            for (int h=1; h<4; ++h) for (int v=0; v<4; ++v)
                if (even[h] & (1<<(4*v+(v^h)))) odd_zero |= 1<<(4*(v^h)+v);
            for (int d=0; d<16; ++d) {
                int diagonal=0;
                for (int v=0; v<4; ++v) if (d & (1<<v)) diagonal |= 1<<(5*v);
                even[0]=even_zero|diagonal;
                odd[0]=odd_zero|diagonal;
                if (!shapes[even[0]].rank[3] || !shapes[odd[0]].rank[3]) continue;
                ++candidates;
                std::array<int,16> edges{};
                for (int p=0; p<4; ++p) for (int h=0; h<4; ++h) {
                    auto& shape=shapes[p%2 ? odd[h] : even[h]];
                    int free=p<2 ? shape.free_rows : shape.free_columns;
                    for (int v=0; v<4; ++v) if (free & (1<<v))
                        edges[4*p+v] |= 1<<(4*((p+1)%4)+(v^h));
                }
                if (has_branching_cycle(edges)) continue;
                ++finite_paths;
                if ((i&j&k&a.mask&b.mask&c.mask)==65535) { ++dense; continue; }
                std::uint64_t mask=0;
                for (int h=0; h<4; ++h) mask |= std::uint64_t(even[h])<<(16*h);
                std::cout<<mask<<'\n';
            }
        }
    }
    if (candidates!=39461472 || finite_paths!=5160 || dense!=16)
        throw std::logic_error("Unexpected (3,1,1,1) counts");
}
