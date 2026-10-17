#pragma once
namespace sdbd::cli {
// Parse run/generate/validate and informational commands from argc/argv; perform requested
// operation and return status, or throw with diagnostics.
int command_line(int argc,char** argv);
}
