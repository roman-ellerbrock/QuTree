#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cerr
        << "Provide input as first argument (and no further arguments).\n";
    exit(1);
  }

  std::string filename(argv[1]);

//  tetrachem::run(filename);

  return EXIT_SUCCESS;
}


