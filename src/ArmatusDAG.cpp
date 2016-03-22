/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
    See LICENSE.txt included with this distribution.
*/

#include <limits>
#include <boost/heap/binomial_heap.hpp>
#include <algorithm>
#include "ArmatusDAG.hpp"


/*

OPT(l) = max{  max_{ k<l } OPTD(k-1),  // OPT(l) ends in a non-domain
               OPTD(l) }               // OPT(l) ends in a domain

OPTD(l) = max_{ k<l } OPT(k-1) + q(k,l)

q(k,l) = {  s(k,l)   if s(k,l) > 0,
            -inf     otherwise  }

OPT(0) = OPT(1) = OPTD(0) = OPTD(1) = 0

*/

ArmatusDAG::ArmatusDAG(ArmatusParams& p) : 
    OPT(SubProbMatrix(p.n+1)),
    OPTD(SubProbMatrix(p.n+1)) { 
	params = &p;
    for (size_t l=0; l <= p.n; l++) {
       OPT[l].resize(p.K); 
       OPTD[l].resize(p.K); 
    }
}

double ArmatusDAG::s(size_t k, size_t l) {
    size_t d_i = d(k,l);
    return params->sums(k-1, l-1)/ std::pow(static_cast<double>(d_i),params->gamma);
}


/*
q(k,l) = {  s(k,l)   if s(k,l) > 0,
            -inf     otherwise }
*/
double ArmatusDAG::q(size_t k, size_t l) {
    size_t d_i = d(k,l);
    double score = (s(k, l) - params->mu[d_i]);
    if (score > 0) {
        return score;
    }
	return -std::numeric_limits<double>::infinity();
}

void ArmatusDAG::build() {
    // OPTD(0) = OPTD(1) = 0 for all K near-optimal solutions
    OPT[0][0] = OPT[1][0] = {0, 1, 0};
    OPTD[0][0] = OPTD[1][0] = {0, 1, 0};
    for (size_t i=1; i<params->K; i++) {
        OPT[0][i] = OPT[1][i] = {-std::numeric_limits<double>::infinity(), 1, 0};
        OPTD[0][i] = OPTD[1][i] = {-std::numeric_limits<double>::infinity(), 1, 0};
    }

    // Build optimal solutions for l=2 to n
    for (size_t l=2; l<=params->n; l++) {

        // Initialize best scores and backpointers
        double scoreDomain, scoreNonDomain;
        size_t backPointerDomain, backPointerNonDomain;

        scoreDomain = scoreNonDomain = -std::numeric_limits<double>::infinity();
        backPointerDomain = backPointerNonDomain = 1;


        // max_{ k<l } OPTD( k-1 )
        for (size_t k=1; k<l; k++) {
            if (OPTD[k-1][0].score > scoreNonDomain) {
                scoreNonDomain = OPTD[k-1][0].score;
                backPointerNonDomain = k;
            }
        }

        // OPTD(l) = max_{ k<l } OPT(k-1) + q(k,l)
        for (size_t k=1; k<l; k++) {
            double candidateScore = OPT[k-1][0].score + q(k,l);
            if (candidateScore > scoreDomain) {
                scoreDomain = candidateScore;
                backPointerDomain = k;
            }
        }
        OPTD[l][0] = {scoreDomain, backPointerDomain, 0};

        /*  OPT(l) = max{  max_{ k<l } OPTD(k-1),  // OPT(l) ends in a non-domain
                           OPTD(l) }               // OPT(l) ends in a domain     */
        if (scoreNonDomain > scoreDomain) {
            OPT[l][0] = {scoreNonDomain, backPointerNonDomain, 0};
        } else {
            OPT[l][0] = OPTD[l][0];
        }
    }
    cout << "OPTIMAL SCORE: " << OPT[params->n][0].score << endl;
}

void ArmatusDAG::computeTopK() {
    std::cerr << "begin computeTopK()\n";
    for (size_t l=2; l<=params->n; l++) {
        /*
        
        OPT(l) = max{  max_{ k<l } OPTD(k-1),  // OPT(l) ends in a non-domain
                       OPTD(l) }               // OPT(l) ends in a domain
        
        OPTD(l) = max_{ k<l } OPT(k-1) + q(k,l)
        
        */
        using heapT = boost::heap::binomial_heap<SubProblem>;
        heapT heapOPTNonDomain;
        heapT heapOPTDomain;

        for (size_t k=1; k<l; k++) {
            auto nonDomainCandidate = OPTD[k-1][0];
            nonDomainCandidate.backPointer = k;
            heapOPTNonDomain.push(nonDomainCandidate); 

            auto domainCandidate = OPT[k-1][0];
            domainCandidate.score += q(k,l);
            domainCandidate.backPointer = k;
            heapOPTDomain.push(domainCandidate);
        }

        auto pushSubProblem = [&](heapT & heap, SubProblem & subProb, bool isDomain) {
            heap.pop();
            auto nextOptimalIndex = subProb.backOptimalIndex + 1;
            auto k = subProb.backPointer;

            if (nextOptimalIndex < params->K) {
                SubProblem newCandidate;
                if (isDomain) {
                    newCandidate = OPT[k-1][nextOptimalIndex];
                    newCandidate.score += q(k,l);
                } else {
                    newCandidate = OPTD[k-1][nextOptimalIndex];
                }
                newCandidate.backPointer = k;
                newCandidate.backOptimalIndex = nextOptimalIndex;
                heap.push(newCandidate);
            }
        };

        size_t i = 0;
        size_t j = 0;
        while (i < params->K) {
            OPTD[l][i] = heapOPTDomain.top();
            pushSubProblem(heapOPTDomain, OPTD[l][i], true);
            auto nonDomainCandidate = heapOPTNonDomain.top();
             
            if (nonDomainCandidate.score > OPTD[l][j].score) {
                OPT[l][i] = nonDomainCandidate;
                pushSubProblem(heapOPTNonDomain, nonDomainCandidate, false);
            } else {
                OPT[l][i] = OPTD[l][j];
                j++;
            }

            i++;    
        }
    }

    std::cerr << "In topK()\n";
    for (size_t i = 0; i < params->K; ++i) {
        std::cerr << "The " << i << "th-best solution had score " << OPT[params->n][i].score << "\n";
    }
}

DomainSet ArmatusDAG::extractDomains(size_t i) {
    size_t k,l;
    DomainSet dset;

    l = params->n;
    do {
        k = OPT[l][i].backPointer;
        if (q(k,l) > 0) {
            dset.push_back(Domain(k-1,l-1)); 
        }
        l = k-1;
    } while(l > 1);

    return dset;
}

WeightedDomainEnsemble ArmatusDAG::extractTopK() {
    WeightedDomainEnsemble ensemble;

    for (size_t i = 0; i < params->K; i++) { 
        ensemble.domainSets.push_back( extractDomains(i) );
        ensemble.weights.push_back( OPT[params->n][i].score/OPT[params->n][0].score );
    }

    return ensemble;
}
