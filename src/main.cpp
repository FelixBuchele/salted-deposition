#include "cli.hpp"
#include <exception>
#include <iostream>
// Run the requested CLI command; print errors and return a nonzero status on failure.
int main(int argc, char** argv) {
    try {
        return sdbd::cli::command_line(argc,argv);
    }
    catch (const std::exception& error) { std::cerr << "deposit: " << error.what() << "\n"; return 1; }
}
