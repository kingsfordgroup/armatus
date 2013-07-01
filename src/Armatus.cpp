#include "Version.hpp"
#include "MatrixParser.hpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

  cout << "***********************************" << endl;
  cout << "***********************************" << endl;
  cout << "****                           ****" << endl;
  cout << "****        ARMATUS 1.0        ****" << endl;
  cout << "****                           ****" << endl;
  cout << "***********************************" << endl;
  cout << "***********************************" << endl;

  MatrixParser parser;
  parser.parseGZipMatrix("SOME_PATH");
  return 0;
}