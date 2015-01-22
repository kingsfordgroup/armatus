/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
*/

#ifndef __ARMATUS_PARAMS_HPP__
#define __ARMATUS_PARAMS_HPP__

#include <memory>

#include <boost/numeric/ublas/symmetric.hpp>

#include "ArmatusUtil.hpp"

class ArmatusParams {
  using SymmetricMatrix = boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>;
  public:
    explicit ArmatusParams(std::shared_ptr<SparseMatrix> _A, double _gamma, size_t Kp, int minMeanSamples);
    void computeSumMuSigma_();
    std::shared_ptr<SparseMatrix> A;
    SymmetricMatrix sums;
    size_t n;
    size_t K;
    size_t minMeanSamples;
    double gamma;
    std::vector<double> mu;
};

#endif // __ARMATUS_PARAMS_HPP__
