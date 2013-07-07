#include <limits>
#include <boost/heap/binomial_heap.hpp>
#include <boost/range/irange.hpp>
#include <algorithm>
#include "ArmatusDAG.hpp"

ArmatusDAG::ArmatusDAG(ArmatusParams& p) :
	subProbs(SubProbsVec(p.n)), 
    edgeWeights(EdgeWeightsVec(p.n, vector<double>())) {
	params = &p;
}

void ArmatusDAG::build() {
	// Since building the DAG visits everything in the same order as computing
	// a single optimial solution, I've had this function also fill out the
	// necessary information to obtain a backtrace of the opt. sln.
	subProbs[0].topK.push_back({std::numeric_limits<size_t>::max(), 0, 0.0});
	// What's the right condition here?
	subProbs[1].topK.push_back({std::numeric_limits<size_t>::max(), 0, q(0,1)});	

	for (size_t l : boost::irange(size_t{2}, params->n)) {
  	    edgeWeights[l].resize(l-1);
  	    size_t chosenEdge = 0;
  	    double bestScore = 0.0;
  	    for (size_t k : boost::irange(size_t{1}, l)) {
  	    	// The last domain is [k,l]
  	    	double edgeWeight = std::max(q(k, l), 0.0);
  	    	// If the last domain starts at k (and is inclusive), then the
  	    	// previous domain ends at k-1, so the backpointer should use the
  	    	// edge to k-1 and the best solution to that subproblem
	    	edgeWeights[l][k-1] = edgeWeight;
	    	// The score of this solution is the benefit we get from traversing
	    	// the edge + the best score available at the child
  	    	double thisScore = edgeWeight + subProbs[k-1].topK[0].score;
  	    	
  	    	// If this is the best derivation of this subproblem so far
  	    	// then record it.
  	    	if (thisScore > bestScore) {
 				bestScore = thisScore;
 				chosenEdge = k-1;
  	    	}

  	    }

  	    // Put our best solution on the stack.
  	    subProbs[l].topK.push_back({chosenEdge, 0, bestScore});
    }
}


double ArmatusDAG::s(size_t k, size_t l) {
	size_t d = l-k;
    return params->sums(k, l)/ std::pow(static_cast<double>(d),params->gamma);
}

double ArmatusDAG::q(size_t k, size_t l) {
	size_t d = l-k;
	return (s(k,l) - params->mu[d]);
}

/**
*  This interface is *wrong* but I wanted to sketch the
*  basic algo.
**/
void ArmatusDAG::computeTopK(uint32_t k) {

	for (size_t l : boost::irange(size_t{2}, params->n)) {
		auto BackPointerComparator = [l, this] (const BackPointer& x, const BackPointer& y) -> bool {
			return (this->edgeWeights[l][x.edge] + this->subProbs[x.edge].topK[x.childSolution].score) < 
			       (this->edgeWeights[l][y.edge] + this->subProbs[y.edge].topK[y.childSolution].score);
		};


		boost::heap::binomial_heap<BackPointer,
		                           boost::heap::compare<decltype(BackPointerComparator)>> pq(BackPointerComparator);
		for (size_t k : boost::irange(size_t{1}, l)) { 
			double score = edgeWeights[l][k-1] + subProbs[k-1].topK[0].score;
			pq.push({k-1, 0, score});
		}

		while (subProbs[l].topK.size() < k and !pq.empty()) {
			auto& bp = pq.top();
			pq.pop();
			if (bp != subProbs[l].topK.back()) { subProbs[l].topK.push_back(bp); }
			size_t nextSlnIdx = bp.childSolution + 1;
			if (nextSlnIdx < subProbs[bp.edge].topK.size()) {
				double score = edgeWeights[l][k-1] + subProbs[k-1].topK[nextSlnIdx].score;
				pq.push({k-1, nextSlnIdx, score});
			}
		}
	}

	size_t sln = 0;
	while (sln < k and sln < subProbs[params->n-1].topK.size()) {
		//std::cerr << "solution " << sln << " has score " << subProbs[params->n-1].topK[sln].score << "\n";
		++sln;
	}
}

vector<Domain> ArmatusDAG::viterbiPath() {
	vector<Domain> domains;
	const size_t INVALID = std::numeric_limits<size_t>::max();
	size_t end = params->n-1;
	//std::cerr << "Best score is " << subProbs[end].topK[0].score << "\n";

	while (subProbs[end].topK[0].edge != INVALID) {
		size_t begin = subProbs[end].topK[0].edge;
		domains.push_back({begin+1, end});
		end = begin;
	}

    sort(domains.begin(), domains.end());

	return domains;
}
