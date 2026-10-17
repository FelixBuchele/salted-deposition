#include "cli.hpp"
#include <exception>
#include <iostream>
// Execute the configuration-first frontend; diagnose handled errors and return nonzero.
int main(int argc, char** argv) {
    try {
        return sdbd::cli::command_line(argc,argv);
    }
    catch (const std::exception& error) { std::cerr << "deposit: " << error.what() << "\n"; return 1; }
}
