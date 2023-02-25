#include <cstddef>

template<typename T> static std::size_t partition(T* a, T* b, T* c, std::size_t* order,
                                                  std::size_t begin, std::size_t end, T mid_value)
{
    std::size_t pivot_idx;
    
    if (begin + 1 < end) {

        T ta = a[begin];
        T tb = b[begin];
        T tc = c[begin];
        
        std::size_t t_idx = order[begin];
        std::size_t upper = begin;
        std::size_t lower = end - 1;
        
        a[begin] = mid_value;
        
        while (upper != lower) {
            while ((upper < lower) && (mid_value < a[lower])) lower--;
            
            if (upper != lower) {
                a[upper] = a[lower];
                b[upper] = b[lower];
                c[upper] = c[lower];
                order[upper] = order[lower];
            }
            
            while ((upper < lower) && (mid_value >= a[upper])) upper++;
            
            if (upper != lower) {
                a[lower] = a[upper];
                b[lower] = b[upper];
                c[lower] = c[upper];
                order[lower] = order[upper];
            }
        }
        
        pivot_idx = upper + 1;
        
        if (ta > mid_value) pivot_idx = upper;
        
        a[upper] = ta;
        b[upper] = tb;
        c[upper] = tc;
        order[upper] = t_idx;
        
    } else if (begin + 1 == end) {
        
        if (a[begin] < mid_value)
            pivot_idx = begin + 1;
        else
            pivot_idx = begin;
        
    } else {
        
        pivot_idx = begin;
    }
    
    return pivot_idx;
}
