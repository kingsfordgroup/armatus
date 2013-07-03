#include "ArmatusParams.hpp"

ArmatusParams::ArmatusParams(std::shared_ptr<SparseMatrix> Ap, double gammap) :
 A(Ap), n(Ap->size1()), gamma(gammap), mu(std::vector<double>(Ap->size1())),
 sigma(std::vector<double>(Ap->size1())), median(std::vector<double>(Ap->size1())),
 sums(ArmatusParams::SymmetricMatrix(Ap->size1(), Ap->size2())) {
 //computeSumMuSigma();
 //computeQmax();
}
