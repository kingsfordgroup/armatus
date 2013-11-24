/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
*/

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
        // was this
  	    // edgeWeights[l].resize(l-1);
        // now this
  	    edgeWeights[l].resize(l);
  	    size_t chosenEdge = 0;
  	    double bestScore = 0.0;
  	    for (size_t k : boost::irange(size_t{0}, l)) {
  	    	// The last domain is [k,l]
  	    	double edgeWeight = std::max(q(k, l), 0.0);
            edgeWeights[l][k] = edgeWeight;

            if (k > 0) {
                // The score of this solution is the benefit we get from traversing
                // the edge + the best score available at the child
                double thisScore = edgeWeight + subProbs[k-1].topK[0].score;

                // If this is the best derivation of this subproblem so far
                // then record it.
                if (thisScore > bestScore) {
                    bestScore = thisScore;
                    chosenEdge = k;
                }
            } else {
                double thisScore = edgeWeight;
                if (thisScore > bestScore) {
                    bestScore = thisScore;
                    chosenEdge = 0;
                }
            }
  	    }

  	    // Put our best solution on the stack.
  	    subProbs[l].topK.push_back({chosenEdge, 0, bestScore});
    }
}


double ArmatusDAG::s(size_t k, size_t l) {
	size_t d = l-k+1;
    return params->sums(k, l)/ std::pow(static_cast<double>(d),params->gamma);
}

double ArmatusDAG::q(size_t k, size_t l) {
	size_t d = l-k+1;
	return (s(k, l) - params->mu[d]);
}

/**
*  This interface is *wrong* but I wanted to sketch the
*  basic algo.
**/

/*
class BackPointerComparator {
    bool operator()(const BackPointer& x, const BackPointer& y) {
	    return (this->edgeWeights[l][x.edge] + this->subProbs[x.edge].topK[x.childSolution].score) < 
    	       (this->edgeWeights[l][y.edge] + this->subProbs[y.edge].topK[y.childSolution].score);
    }
};
*/

void ArmatusDAG::computeTopK(uint32_t k) {
	for (size_t l : boost::irange(size_t{2}, params->n)) {
		std::function<bool (const BackPointer&, const BackPointer&)> BackPointerComparator = [l, this] (const BackPointer& x, const BackPointer& y) -> bool {
            auto scoreX = this->edgeWeights[l][x.edge];
            if (x.edge > 0) { scoreX += this->subProbs[x.edge-1].topK[x.childSolution].score; }
            auto scoreY = this->edgeWeights[l][y.edge];
            if (y.edge > 0) { scoreY += this->subProbs[y.edge-1].topK[y.childSolution].score; }

			return scoreX < scoreY;
		};

		boost::heap::binomial_heap<BackPointer,
		                           boost::heap::compare<decltype(BackPointerComparator)>> pq(BackPointerComparator);
		for (size_t j : boost::irange(size_t{0}, l)) { 
            if (j > 0) {
                size_t i = j - 1;
                bool isDomain = edgeWeights[l][j] > 0; 
                bool prevIsDomain = (subProbs[i].topK.size() > 0 and edgeWeights[i].size() > 0) ? edgeWeights[i][subProbs[i].topK[0].edge] > 0 : true;

                if (isDomain or prevIsDomain) {
                    double score = edgeWeights[l][j] + subProbs[i].topK[0].score;
                    pq.push({j, 0, score});
                }
            } else {
                double score = edgeWeights[l][0];
                pq.push({0, 0, score});
            }
		}

		while (subProbs[l].topK.size() < k and !pq.empty()) {
			auto bp = pq.top();
			pq.pop();
            size_t j = bp.edge;
            if (j > 0) {
                if (bp != subProbs[l].topK.back()) { subProbs[l].topK.push_back(bp); }
                size_t nextSlnIdx = bp.childSolution + 1;
                bool isDomain = edgeWeights[l][j] > 0;
                
                size_t i = j - 1;
                //std::cout << bp.edge << "\t" << subProbs.size() << endl;
                if (nextSlnIdx < subProbs[i].topK.size()) {
                    bool isPrevDomain = edgeWeights[i][subProbs[i].topK[nextSlnIdx].edge] > 0;
                    if (isDomain or isPrevDomain) {
                        double score = edgeWeights[l][j] + subProbs[i].topK[nextSlnIdx].score;
                        pq.push({j, nextSlnIdx, score});
                    }
                }
            } else {
                if (bp != subProbs[l].topK.back()) { subProbs[l].topK.push_back(bp); }
            }
		}
	}

	// size_t sln = 0;
	// while (sln < k and sln < subProbs[params->n-1].topK.size()) {
	// 	std::cerr << "solution " << sln << " has score " << subProbs[params->n-1].topK[sln].score << "\n";
	// 	++sln;
	// }
}


WeightedDomainEnsemble ArmatusDAG::extractTopK(uint32_t k) {

  const size_t INVALID = std::numeric_limits<size_t>::max();
  // We skip 0 since that was already extracted by viterbiPath
  uint32_t currSln{0};

  WeightedDomainEnsemble solutions{ DomainEnsemble(k, DomainSet()), Weights(k, 0.0) };
  auto root = params->n-1;
  auto topScore = subProbs[root].topK[0].score;
  auto scores = Weights(k, 0.0);

  while (currSln < k and currSln < subProbs[root].topK.size()) {
    auto currentScore = subProbs[root].topK[currSln].score;
    scores[currSln] = currentScore;
    auto& currentDomainSet = solutions.domainSets[currSln];
    solutions.weights[currSln] = currentScore / topScore;

    // Which solution to use at the child
    auto bp = currSln;
    int64_t end = params->n-1;

    while (end >= 0 and subProbs[end].topK[bp].edge != INVALID) {
      size_t begin = subProbs[end].topK[bp].edge;
      size_t start = begin;
//      if (subProbs[begin].topK[subProbs[end].topK[bp].childSolution].edge != INVALID)
//         start++;
      Domain d(start, end);
      if (d.score(*params) > 0.0){
         currentDomainSet.push_back(d);
      }
      bp = subProbs[end].topK[bp].childSolution;
      end = begin-1;
    }

    ++currSln;
  }

  /**
   * consistency check all the extracted domains (for debug purposes only, 
   * take this out for release / speed).
   */
  size_t c = 0;
  double eps = 1e-3;
  for (auto& ds : solutions.domainSets) {
    double score = 0.0;
    for (auto& d : ds) {
      score += d.score(*params);
    }
    if (std::abs(score - scores[c]) > eps) { 
      std::cerr << "For solution [" << c << "], the recorded score was " << scores[c] << " but extracted domains sum to " << score << "\n";
    }
    ++c;
  }
  
  return solutions;
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
