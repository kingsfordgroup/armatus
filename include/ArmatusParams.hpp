#ifndef __ARMATUS_PARAMS_HPP__
#define __ARMATUS_PARAMS_HPP__

#include <memory>

#include <boost/numeric/ublas/symmetric.hpp>

#include "MatrixParser.hpp"

class ArmatusParams {
  using SymmetricMatrix = boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>;
  public:
    explicit ArmatusParams(std::shared_ptr<MatrixParser::SparseMatrix> _A, double _gamma);
    void computeSumMuSigma_();
    std::shared_ptr<MatrixParser::SparseMatrix> A;
    SymmetricMatrix sums;
    size_t n;
    double gamma;
    std::vector<double> mu;
};

#endif // __ARMATUS_PARAMS_HPP__