#include <stddef.h>
#include <stdlib.h>
#include <new>

import preprocess;

static size_t allocations;
void* operator new(size_t n) {
    ++allocations;
    if (void* p = malloc(n)) return p;
    throw std::bad_alloc();
}
void* operator new[](size_t n) {
    ++allocations;
    if (void* p = malloc(n)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { free(p); }
void operator delete[](void* p) noexcept { free(p); }
void operator delete(void* p, size_t) noexcept { free(p); }
void operator delete[](void* p, size_t) noexcept { free(p); }

struct foreign_cell { float expression; int label; double custom_feature; };
struct fusion {
    const foreign_cell* in;
    float* out;
    float gain, bias;
    size_t calls;
    void operator()(size_t first, size_t count) {
        ++calls;
        for (size_t i = first; i < first+count; ++i)
            out[i] = gain*in[i].expression + bias*float(in[i].custom_feature);
    }
};
void runtime_fusion(void* raw, size_t first, size_t count) {
    auto& f = *static_cast<fusion*>(raw);
    ++f.calls;
    for (size_t i = first; i < first+count; ++i)
        f.out[i] = f.gain*f.in[i].expression + f.bias*float(f.in[i].custom_feature);
}

int main() {
    const foreign_cell input[] = {{1,8,2},{3,9,4},{5,10,6},{7,11,8},{9,12,10}};
    float a[5] = {}, b[5] = {};
    fusion stat = {input,a,2,.5f,0}, dyn = {input,b,2,.5f,0};
    size_t before = allocations;
    auto s = preprocess::for_blocks(5,2,stat);
    auto r = preprocess::for_blocks(5,2,runtime_fusion,&dyn);
    size_t after = allocations;
    if (!s || !r || before != after) return 1;
    if (stat.calls != 3 || dyn.calls != 3) return 2;
    for (size_t i = 0; i < 5; ++i) {
        float expected = 2*input[i].expression + .5f*float(input[i].custom_feature);
        if (a[i] != expected || b[i] != expected) return 3;
    }
    if (preprocess::for_blocks(5,0,stat).code != preprocess::QcCode::invalid_argument) return 4;
    if (preprocess::for_blocks(5,2,nullptr,&dyn).code != preprocess::QcCode::invalid_argument) return 5;
    auto empty = preprocess::for_blocks(0,2,stat);
    if (!empty || stat.calls != 3) return 6;
    return 0;
}
