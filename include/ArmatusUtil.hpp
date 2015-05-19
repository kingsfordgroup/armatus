/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
*/

#ifndef __ARMATUS_UTIL_HPP__
#define __ARMATUS_UTIL_HPP__

#include <memory>
#include <string>
#include <vector>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;

//using SparseMatrix = boost::numeric::ublas::compressed_matrix<double>;
using SparseMatrix = boost::numeric::ublas::matrix<double>;

class MatrixProperties {
    public:
    std::shared_ptr<SparseMatrix> matrix;
    std::string chrom;
    int resolution;
};

MatrixProperties parseGZipMatrix(string path, int resolution, string chrom);
MatrixProperties parseRaoMatrix(string path, int resolution, string chrom, bool noNormalization);

double d(size_t const & i, size_t const & j);

class ArmatusParams;
class Domain {
    public:
    Domain(size_t s, size_t e);
    size_t start;
    size_t end;
    double score(ArmatusParams& p);
    bool operator< (const Domain &other) const {
        return std::make_tuple(start, end) < std::make_tuple(other.start, other.end);
    }
};

using DomainSet = vector<Domain>;
using DomainEnsemble = vector<DomainSet>;
using Weights = vector<double>;
using Resolutions = vector<double>;
using OptIndex = vector<size_t>;

struct WeightedDomainEnsemble {
    DomainEnsemble domainSets;
    Weights weights;
    Resolutions resolutions;
    OptIndex optidx;
};

// void optimalDomains(std::shared_ptr<SparseMatrix> A, float gamma);

WeightedDomainEnsemble multiscaleDomains(std::shared_ptr<SparseMatrix> A, float gammaMax, double stepSize, int k, int minMeanSamples, bool justThisGamma);
DomainSet consensusDomains(WeightedDomainEnsemble& dEnsemble);

void outputDomains(DomainSet dSet, string fname, MatrixProperties matProp);


#endif // __ARMATUS_UTIL_HPP__
