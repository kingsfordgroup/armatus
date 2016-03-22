/* 
	Authors: Darya Filippova, Geet Duggal, Rob Patro
	dfilippo | geet | robp @cs.cmu.edu
        See LICENSE.txt included with this distribution.
*/


#include <iomanip>
#include <limits>

#include <fstream>

#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>

#include "ArmatusParams.hpp"

ArmatusParams::ArmatusParams(std::shared_ptr<SparseMatrix> Ap, double gammap, size_t Kp, int minMeanSamples) :
    A(Ap), 
    sums(ArmatusParams::SymmetricMatrix(Ap->size1(), Ap->size2())), 
    n(Ap->size1()), 
    K(Kp), 
    minMeanSamples(minMeanSamples),
    gamma(gammap), 
    mu(std::vector<double>(Ap->size1()+1))
{
 computeSumMuSigma_();
}

void ArmatusParams::computeSumMuSigma_() {
	using namespace boost::accumulators;
	using namespace boost::range;
	// Vector to hold accumulators that will compute the mean
	// for domains of each size.
	std::vector<accumulator_set<double,stats<tag::mean, tag::count>>> acc(n+1);
	//std::vector<accumulator_set<double,stats<tag::median(with_p_square_quantile), tag::count>>> acc(n+1);

	// A reference will be easier to work with here
	SparseMatrix& M = *A;

	for (size_t i : boost::irange(size_t{0}, n)) {
		sums(i, i) = M(i, i);
	}

	for (size_t i : boost::irange(size_t{1}, n)) {
		std::vector<double> columnSums(i+1);
		columnSums[i] = M(i, i);
		for (size_t j : boost::adaptors::reverse(boost::irange(size_t{0}, i))) {
			columnSums[j] = columnSums[j+1] + M(j, i);
			sums(j, i) = sums(j, i-1) + columnSums[j];
			int d_i = d(j,i);
			assert(d_i >= 0);
			double s = sums(j, i) / std::pow(static_cast<double>(d_i), gamma);
		    acc[d_i](s);
		}
	}

	for (size_t i : boost::irange(size_t{0}, n+1)) {
		mu[i] = mean(acc[i]);
        //mu[i] = median(acc[i]);
		// Require at least 100 samples to compute a Z-score
		if (boost::accumulators::count(acc[i]) < minMeanSamples) { 
			mu[i] = std::numeric_limits<double>::max();
		}
	}
}
