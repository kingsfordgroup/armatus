#ifndef __ARMATUS_UTIL_HPP__
#define __ARMATUS_UTIL_HPP__

#include <memory>
#include <string>
#include <vector>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using SparseMatrix = boost::numeric::ublas::compressed_matrix<double>;

std::shared_ptr<SparseMatrix> parseGZipMatrix(std::string path);

class ArmatusParams;
class Domain {
    public:
    Domain(size_t s, size_t e);
    size_t start;
    size_t end;
    double score(ArmatusParams& p);
    bool operator< (const Domain &other) const {
        return end < other.end;
    }
};

using DomainSet = std::vector<Domain>;
using DomainEnsemble = std::vector<DomainSet>;

DomainSet consensus(DomainEnsemble dEnsemble);


#endif // __ARMATUS_UTIL_HPP__