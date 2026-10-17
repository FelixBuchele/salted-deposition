#pragma once
#include <iostream>
namespace sdbd::detail {
// Print diagnostic values separated by three spaces, followed by a flushed newline.
// Keep the original formatting for runs using positional arguments.
template<class First, class... Rest>
void ds(const First& first, const Rest&... rest) {
    std::cout << first;
    ((std::cout << "   " << rest), ...);
    std::cout << std::endl;
}
}
