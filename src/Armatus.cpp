#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "ArmatusParams.hpp"
#include "Version.hpp"
#include "ArmatusUtil.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  auto str = R"(
	***********************************
	***********************************
	****                           ****
	****        ARMATUS 1.0        ****
	****                           ****
	***********************************
	***********************************
	)";

  namespace po = boost::program_options;
  using std::string;
  using std::cerr;

  // Declare the supported options.
  po::options_description opts("armatus options");
  opts.add_options()
  ("help,h", "produce help message")
  ("input,i", po::value<string>()->required(), "input file")
  ("gamma,g", po::value<double>()->required(), "gamma [scaling] parameter")
  ;

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opts), vm);

    if (vm.count("help")) {
      cerr << str << "\n";
      cerr << opts << "\n";
      std::exit(1);
    }

    po::notify(vm);    

    if (vm.count("input")) {
      cerr << "Reading input from " << vm["input"].as<string>() << ".\n";
      // parse input matrix file

  	  auto mat = parseGZipMatrix(vm["input"].as<string>());
      cerr << "MatrixParser read matrix of size: " << mat->size1() << " x " << mat->size2()  << "\n";
      ArmatusParams params(mat, vm["gamma"].as<double>());
    } else {
      cerr << "Input file was not set.\n";
      std::exit(1);
    }

  } catch (po::error& e) {
    cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  }

  

  return 0;
}