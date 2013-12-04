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
    OPT(vector<SubProblem>(p.n+1)),
    OPTD(vector<SubProblem>(p.n+1)) { 
	params = &p;
    for (size_t l=0; l <= p.n; l++) {
       OPT[l].score.resize(p.K); 
       OPTD[l].score.resize(p.K); 
    }
}

double ArmatusDAG::s(size_t k, size_t l) {
	size_t d = l-k+1;
    return params->sums(k-1, l-1)/ std::pow(static_cast<double>(d),params->gamma);
}


/*
q(k,l) = {  s(k,l)   if s(k,l) > 0,
            -inf     otherwise }
*/
double ArmatusDAG::q(size_t k, size_t l) {
	size_t d = l-k+1;
    double score = (s(k, l) - params->mu[d]);
    if (score > 0) return score;
	return -std::numeric_limits<double>::infinity();
}

void ArmatusDAG::build() {
    // OPTD(0) = OPTD(1) = 0 for all K near-optimal solutions
    for (size_t i=0; i<params->K; i++) {
        OPT[0].score[i] = OPT[1].score[i] = 0;
        OPTD[0].score[i] = OPTD[1].score[i] = 0;
    }

    // Build optimal solutions for l=2 to n
    for (size_t l=2; l<=params->n; l++) {
        double scoreDomain, scoreNonDomain;
        scoreDomain = scoreNonDomain = -std::numeric_limits<double>::infinity();

        // max_{ k<l } OPTD( k-1 )
        for (size_t k=1; k<l; k++) {
            if (OPTD[k-1].score[0] > scoreNonDomain) {
                scoreNonDomain = OPTD[k-1].score[0];
            }
        }

        // OPTD(l) = max_{ k<l } OPT(k-1) + q(k,l)
        for (size_t k=1; k<l; k++) {
            //cout << "k=" << k << " l=" << l << " q(k,l)=" << q(k,l) << endl;
            double candidateScore = OPT[k-1].score[0] + q(k,l);
            if (candidateScore > scoreDomain) {
                scoreDomain = candidateScore;
            }
        }
        OPTD[l].score[0] = scoreDomain;

        /*  OPT(l) = max{  max_{ k<l } OPTD(k-1),  // OPT(l) ends in a non-domain
                           OPTD(l) }               // OPT(l) ends in a domain     */
        if (scoreNonDomain > scoreDomain) {
            OPT[l].score[0] = scoreNonDomain;
        } else {
            OPT[l].score[0] = scoreDomain;
        }
    }
    cout << "OPTIMAL SCORE: " << OPT[params->n].score[0] << endl;
}

void ArmatusDAG::computeTopK(uint32_t k) {
}


WeightedDomainEnsemble ArmatusDAG::extractTopK(uint32_t k) {
}

