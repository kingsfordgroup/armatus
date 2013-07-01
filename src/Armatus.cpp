#include <boost/program_options.hpp>
#include "Version.hpp"

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;
  using std::string;
  using std::cerr;

  // Declare the supported options.
  po::options_description opts("armatus options");
  opts.add_options()
  ("help", "produce help message")
  ("i,input", po::value<string>(), "input file")
  ;

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);    

    if (vm.count("help")) {
      cerr << opts << "\n";
      return 1;
    }

    if (vm.count("input")) {
      cerr << "Reading input from" << vm["input"].as<string>() << ".\n";
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