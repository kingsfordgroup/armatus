#include <iomanip>
#include <limits>

#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "ArmatusParams.hpp"

ArmatusParams::ArmatusParams(std::shared_ptr<SparseMatrix> Ap, double gammap) :
 A(Ap), n(Ap->size1()), gamma(gammap), mu(std::vector<double>(Ap->size1())),
 sums(ArmatusParams::SymmetricMatrix(Ap->size1(), Ap->size2())) {
 computeSumMuSigma_();
}

void ArmatusParams::computeSumMuSigma_() {
	using namespace boost::accumulators;
	using namespace boost::range;
	// Vector to hold accumulators that will compute the mean
	// for domains of each size.
	std::vector<accumulator_set<double,stats<tag::immediate_mean, tag::count>>> acc(n);

	// A reference will be easier to work with here
	SparseMatrix& M = *A;
	std::cerr << "M is " << M.size1() << " x " << M.size2() << "\n";

	for (size_t i : boost::irange(size_t{0}, n)) {
		sums(i, i) = M(i, i);
	}

	for (size_t i : boost::irange(size_t{1}, n)) {
		std::vector<double> columnSums(i+1);
		columnSums[i] = M(i, i);
		for (size_t j : boost::adaptors::reverse(boost::irange(size_t{0}, i))) {
			columnSums[j] = columnSums[j+1] + M(j, i);
			sums(j, i) = sums(j, i-1) + columnSums[j];
			//sums(i, j) = sums(j, i);

			int d = i - j;
			double s = sums(j, i) / std::pow(static_cast<double>(d), gamma);
			acc[d](s);
		}
	}
		
	for (size_t i : boost::irange(size_t{0}, n)) {
		mu[i] = mean(acc[i]);
		// Require at least 100 samples to compute a Z-score
		if (count(acc[i]) < 100) { 
			mu[i] = std::numeric_limits<double>::max();
		}
	}

}
