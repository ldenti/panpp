#include <cstring>
#include <iostream>

using namespace std;

int main_index(int argc, char *argv[]);
int main_search(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr << "XXX1" << endl;
    // cerr << MAIN_USAGE << endl;
    return 1;
  }
  if (strcmp(argv[1], "index") == 0)
    return main_index(argc - 1, argv + 1);
  else if (strcmp(argv[1], "search") == 0)
    return main_search(argc - 1, argv + 1);

  cerr << "XXX2" << endl;
  // cerr << MAIN_USAGE << endl;
  return 1;
}
