/* 
  Authors: Darya Filippova, Geet Duggal, Rob Patro
  dfilippo | geet | robp @cs.cmu.edu
*/


#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "Version.hpp"
#include "ArmatusUtil.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  auto str = R"(
	***********************************
	****                           ****
	****        ARMATUS 1.0        ****
	****                           ****
	***********************************
	)";

  namespace po = boost::program_options;
  using std::string;
  using std::cerr;


  class Params {
    public:
    string inputFile;
    string outputPrefix;
    double gammaMax;
    size_t k;
    double stepSize;
    bool outputMultiscale;
  };

  Params p;

  // Declare the supported options.
  po::options_description opts("armatus options");
  opts.add_options()
  ("help,h", "produce help message")
  ("input,i", po::value<string>(&p.inputFile)->required(), "input file")
  ("gammaMax,g", po::value<double>(&p.gammaMax)->required(), "gamma-max (highest resolution to generate domains)")
  ("output,o", po::value<string>(&p.outputPrefix)->required(), "output filename prefix")
  ("topK,k", po::value<size_t>(&p.k)->default_value(1), "Compute the top k optimal solutions")
  ("stepSize,s", po::value<double>(&p.stepSize)->default_value(0.05), "Step size to increment resolution parameter")
  ("outputMultiscale,m", po::value<bool>(&p.outputMultiscale)->zero_tokens()->default_value(false), "Output multiscale domains to files as well")
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
      if (p.outputMultiscale) cerr << "Multiresoultion ensemble will be written to files" << endl;
      cerr << "Reading input from " << p.inputFile << ".\n";

  	  auto matProp = parseGZipMatrix(p.inputFile);
      auto mat = matProp.matrix;
      cerr << "MatrixParser read matrix of size: " << mat->size1() << " x " << mat->size2()  << "\n";

      auto dEnsemble = multiscaleDomains(mat, p.gammaMax, p.stepSize, p.k);
      auto dConsensus = consensusDomains(dEnsemble);


      auto consensusFile = p.outputPrefix + ".consensus.txt";
      cerr << "Writing consensus domains to: " << consensusFile << endl;;
      outputDomains(dConsensus, consensusFile, matProp);

      if (p.outputMultiscale) {
        double gamma = 0;
        size_t multiOptIdx = 0;
       cerr << "Writing multiscale domains" << endl;
        for (auto dSetIdx : boost::irange(size_t{0}, dEnsemble.domainSets.size())) {
          auto& dSet = dEnsemble.domainSets[dSetIdx];
          stringstream multiscaleFile;
          multiscaleFile << p.outputPrefix << ".gamma." << std::fixed << std::setprecision(log10(p.gammaMax/p.stepSize)+1) << gamma << "." << setfill('0') << setw(log10(p.k)+1) << multiOptIdx << ".txt";
          outputDomains(dSet, multiscaleFile.str(),matProp);
          multiOptIdx++;
          if (multiOptIdx % p.k == 0)  {
            gamma += p.stepSize;
            multiOptIdx = 0;
          }
        }
      }


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
