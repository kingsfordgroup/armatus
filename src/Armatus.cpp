/* 
  Authors: Darya Filippova, Geet Duggal, Rob Patro
  dfilippo | geet | robp @cs.cmu.edu
  See LICENSE.txt included with this distribution.
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
	****        ARMATUS 2.1        ****
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
    int minMeanSamples;
    int resolution;
    bool outputMultiscale;
    bool raoFormat;
    bool noNormalization;  // No normalization, just raw counts (Rao format)
    bool justGammaMax;
    string chrom;
  };

  Params p;

  // Declare the supported options.
  po::options_description opts("armatus options");
  opts.add_options()
  ("parseRaoFormat,R", po::value<bool>(&p.raoFormat)->zero_tokens()->default_value(false), "Parse the Rao data format instead of Dixon et al. (Provide, for example, GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb,  and KR normalization is used.")
  ("noNormalization,N", po::value<bool>(&p.noNormalization)->zero_tokens()->default_value(false), "No normalization (Rao et al. format)")
  ("gammaMax,g", po::value<double>(&p.gammaMax)->required(), "gamma-max (highest resolution to generate domains)")
  ("justGammaMax,j", po::value<bool>(&p.justGammaMax)->zero_tokens()->default_value(false), "Just obtain domains at the maximum Gamma")
  ("help,h", "produce help message")
  ("input,i", po::value<string>(&p.inputFile)->required(), "input matrix file")
  ("topK,k", po::value<size_t>(&p.k)->default_value(1), "Compute the top k optimal solutions")  
  ("outputMultiscale,m", po::value<bool>(&p.outputMultiscale)->zero_tokens()->default_value(false), "Output multiscale domains to files as well")
  ("resolution,r", po::value<int>(&p.resolution)->default_value(1), "Resolution of data")
  ("chromosome,c", po::value<string>(&p.chrom)->default_value("N/A"), "Chromosome")
  ("minMeanSamples,n", po::value<int>(&p.minMeanSamples)->default_value(100), "Minimum required number of samples to compute a mean")
  ("output,o", po::value<string>(&p.outputPrefix)->required(), "output filename prefix")  
  ("stepSize,s", po::value<double>(&p.stepSize)->default_value(0.05), "Step size to increment resolution parameter")
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

      MatrixProperties matProp;
      if (p.raoFormat) {
          matProp = parseRaoMatrix(p.inputFile, p.resolution, p.chrom, p.noNormalization);
      }
      else { 
          matProp = parseGZipMatrix(p.inputFile, p.resolution, p.chrom);
      }
      auto mat = matProp.matrix;
      cerr << "MatrixParser read matrix of size: " << mat->size1() << " x " << mat->size2()  << "\n";

      auto dEnsemble = multiscaleDomains(mat, p.gammaMax, p.stepSize, p.k, p.minMeanSamples, p.justGammaMax);

      cerr << "Domain ensemble size: " << dEnsemble.domainSets.size() << endl;
      auto dConsensus = consensusDomains(dEnsemble);


      auto consensusFile = p.outputPrefix + ".consensus.txt";
      cerr << "Writing consensus domains to: " << consensusFile << endl;;
      outputDomains(dConsensus, consensusFile, matProp);

      if (p.outputMultiscale) {
        cerr << "Writing multiscale domains" << endl;
        for (auto dSetIdx : boost::irange(size_t{0}, dEnsemble.domainSets.size())) {
          auto& dSet = dEnsemble.domainSets[dSetIdx];
          size_t multiOptIdx = dEnsemble.optidx[dSetIdx];
          float gamma = dEnsemble.resolutions[dSetIdx];
          stringstream multiscaleFile;
          //multiscaleFile << p.outputPrefix << ".gamma." << std::fixed << std::setprecision(log10(p.gammaMax/p.stepSize)+1) << gamma << "." << setfill('0') << setw(log10(p.k)+1) << multiOptIdx << ".txt";
          multiscaleFile << p.outputPrefix << ".gamma." << gamma << "." << multiOptIdx << ".txt";
          cout << dSetIdx << "\t" << multiscaleFile.str() << endl;
          outputDomains(dSet, multiscaleFile.str(),matProp);
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
