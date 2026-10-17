#pragma once
#include <iostream>
namespace sdbd::detail {
// Print values to stdout with original separators and a flushed newline;
// diagnostics preserve positional-interface output bytes.
template<class First, class... Rest>
void ds(const First& first, const Rest&... rest) {
    std::cout << first;
    ((std::cout << "   " << rest), ...);
    std::cout << std::endl;
}
}
