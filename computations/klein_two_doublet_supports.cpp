// Sparse support orbits for two rank-two and two rank-one colors on each cut.
#include "klein_support_shape.hpp"
#include <limits>
#include <unordered_map>
#include <unordered_set>

using WordMask = std::uint64_t;
using Cuts = std::array<Four,2>;

int bit_index(const Four& word) {
    return 16*(word[0]^word[1])+4*word[0]+word[2];
}

Four word_at(int index) {
    int b=index/16, i=index/4%4, k=index%4;
    return {i,i^b,k,k^b};
}

WordMask block_words(int p, int b, int matrix) {
    WordMask mask=0;
    for (int i=0; i<4; ++i) for (int k=0; k<4; ++k) {
        if (!(matrix & (1<<(4*i+k)))) continue;
        Four word{i,i^b,k,k^b};
        if (p) std::rotate(word.begin(),word.end()-1,word.end());
        mask |= WordMask(1)<<bit_index(word);
    }
    return mask;
}

Cuts slices(WordMask mask) {
    Cuts m{};
    for (WordMask s=mask; s; s&=s-1) {
        auto [i,j,k,l]=word_at(__builtin_ctzll(s));
        m[0][i^j] |= 1<<(4*i+k);
        m[1][j^k] |= 1<<(4*j+l);
    }
    return m;
}

struct Permutation {
    std::array<std::array<WordMask,256>,8> table{};
    WordMask operator()(WordMask mask) const {
        WordMask value=0;
        for (const auto& byte:table) {
            value |= byte[mask&255];
            mask >>= 8;
        }
        return value;
    }
};

std::vector<Permutation> group(int q) {
    std::vector<Permutation> perms;
    for (int rev=0; rev<2; ++rev) for (int rot=0; rot<4; ++rot)
    for (int lin=0; lin<(q==1 ? 2 : 1); ++lin)
    for (int t=0; t<4; ++t) for (int h:{0,1}) for (int k:{0,q}) {
        Permutation pm;
        for (int j=0; j<64; ++j) {
            Four seq=word_at(j);
            if (rev) std::reverse(seq.begin(),seq.end());
            std::rotate(seq.begin(),seq.begin()+rot,seq.end());
            for (int& v:seq) {
                if (q==2 && rot%2) v=Four{0,2,1,3}[v];
                if (lin) v=Four{0,1,3,2}[v];
            }
            seq[0]^=t; seq[1]^=t^h; seq[2]^=t^h^k; seq[3]^=t^k;
            pm.table[j/8][1<<(j%8)]=WordMask(1)<<bit_index(seq);
        }
        for (auto& byte:pm.table) for (int j=1; j<256; ++j)
            byte[j]=byte[j&-j]|byte[j&(j-1)];
        perms.push_back(pm);
    }
    return perms;
}

bool exponential(const Cuts& m, const std::vector<Shape>& shapes) {
    std::array<int,16> edges{};
    for (int p=0; p<4; ++p) for (int b=0; b<4; ++b) {
        const auto& shape=shapes[m[p%2][b]];
        int free=p<2 ? shape.free_rows : shape.free_columns;
        for (int v=0; v<4; ++v) if (free & (1<<v))
            edges[4*p+v] |= 1<<(4*((p+1)%4)+(v^b));
    }
    return has_branching_cycle(edges);
}

int main(int argc, char** argv) {
    if (argc!=2) return 2;
    int q=std::stoi(argv[1]);
    if (q!=1 && q!=2) return 2;
    std::vector<Shape> shapes(65536);
    std::vector<int> rects;
    for (int m=0; m<65536; ++m) {
        shapes[m]=classify(m);
        if (shapes[m].rank[1]) rects.push_back(m);
    }
    WordMask even=block_words(0,2,65535)|block_words(0,3,65535);
    WordMask odd=block_words(1,3^q,65535)|block_words(1,3,65535);
    std::unordered_map<WordMask,std::vector<WordMask>> table;
    for (int i:rects) for (int j:rects) {
        auto w=block_words(1,3^q,i)|block_words(1,3,j);
        table[w&even].push_back(w);
    }
    std::unordered_set<WordMask> thin;
    for (int i:rects) for (int j:rects) {
        auto w=block_words(0,2,i)|block_words(0,3,j);
        for (auto z:table[w&odd]) thin.insert(w|z);
    }
    auto perms=group(q);
    WordMask holes=~(even|odd);
    std::vector<WordMask> cells;
    for (auto s=holes; s; s&=s-1) cells.push_back(s&-s);
    if (cells.size()!=16) throw std::logic_error("Expected 16 missing cells");
    std::vector<WordMask> supplement(65536);
    std::vector<Four> data(65536);
    for (int j=1; j<65536; ++j) {
        supplement[j]=supplement[j&(j-1)]|cells[__builtin_ctz(j)];
        auto m=slices(supplement[j]);
        data[j]={m[0][0],m[0][1],m[1][0],m[1][q]};
    }
    std::uint64_t count=0, kept=0;
    int initial=thin.size(), thin_orbits=0, final_orbits=0;
    while (!thin.empty()) {
        WordMask seed=*thin.begin();
        int weight=0;
        std::vector<int> stabilizer;
        for (int h=0; h<int(perms.size()); ++h) {
            auto z=perms[h](seed);
            weight+=thin.erase(z);
            if (z==seed) stabilizer.push_back(h);
        }
        ++thin_orbits;
        if (seed==(even|odd)) continue;
        auto base=slices(seed);
        for (int j=0; j<65536; ++j) {
            auto m=base;
            m[0][0]|=data[j][0]; m[0][1]|=data[j][1];
            m[1][0]|=data[j][2]; m[1][q]|=data[j][3];
            if (!shapes[m[0][0]].rank[2] || !shapes[m[0][1]].rank[2] ||
                !shapes[m[1][0]].rank[2] || !shapes[m[1][q]].rank[2]) continue;
            count+=weight;
            if (exponential(m,shapes)) continue;
            kept+=weight;
            WordMask w=seed|supplement[j];
            if (std::any_of(stabilizer.begin(),stabilizer.end(),[&](int h){return perms[h](w)<w;})) continue;
            ++final_orbits;
            WordMask canonical=w;
            for (auto& perm:perms) canonical=std::min(canonical,perm(w));
            std::cout<<canonical<<'\n';
        }
    }
    std::array<std::uint64_t,5> got{std::uint64_t(initial),std::uint64_t(thin_orbits),count,kept,
                                    std::uint64_t(final_orbits)};
    std::array<std::uint64_t,5> parallel{1033395,4856,922893330,476008,2779};
    std::array<std::uint64_t,5> transverse{973921,8005,828689160,292056,2794};
    if (got!=(q==1 ? parallel : transverse)) throw std::logic_error("Unexpected support counts");
}
