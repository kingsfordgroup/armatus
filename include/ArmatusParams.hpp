#ifndef __ARMATUS_PARAMS_HPP__
#define __ARMATUS_PARAMS_HPP__

#include <memory>

#include <boost/numeric/ublas/symmetric.hpp>

#include "MatrixParser.hpp"

class ArmatusParams {
  using SymmetricMatrix = boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>;
  public:
    explicit ArmatusParams(std::shared_ptr<MatrixParser::SparseMatrix> _A, double _gamma);

    std::shared_ptr<MatrixParser::SparseMatrix> A;
    SymmetricMatrix sums;
    std::vector<double> mu;
    std::vector<double> sigma;
    std::vector<double> median;
    std::vector<uint32_t> numSamps;
    double qmax;
    double gamma;
    size_t n;
};

#endif // __ARMATUS_PARAMS_HPP__