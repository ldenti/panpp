#include <cstring>
#include <iostream>

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include "usage.hpp"

using namespace std;

int main_index(int argc, char *argv[]);
int main_search(int argc, char *argv[]);
int main_exact(int argc, char *argv[]);
int main_pingpong(int argc, char *argv[]);

int main_fmdindex(int argc, char *argv[]);
int main_fmdexact(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  spdlog::set_default_logger(spdlog::stderr_color_st("stderr"));
  if (argc < 2) {
    cerr << MAIN_USAGE << endl;
    return 1;
  }
  if (strcmp(argv[1], "index") == 0)
    return main_index(argc - 1, argv + 1);
  else if (strcmp(argv[1], "search") == 0)
    return main_pingpong(argc - 1, argv + 1);
  else if (strcmp(argv[1], "exact") == 0)
    return main_exact(argc - 1, argv + 1);
  else if (strcmp(argv[1], "test") == 0)
    return main_search(argc - 1, argv + 1);
  else if (strcmp(argv[1], "fmdindex") == 0)
    return main_fmdindex(argc - 1, argv + 1);
  else if (strcmp(argv[1], "fmdexact") == 0)
    return main_fmdexact(argc - 1, argv + 1);
  else if (strcmp(argv[1], "version") == 0) {
    cerr << VERSION << endl;
    return 0;
  }
  cerr << MAIN_USAGE << endl;
  return 1;
}
