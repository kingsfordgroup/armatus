/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
*/

#ifndef __ARMATUS_DAG_HPP__
#define __ARMATUS_DAG_HPP__

#include "ArmatusUtil.hpp"
#include "ArmatusParams.hpp"
#include <vector>

using namespace boost;
using namespace std;

class ArmatusDAG {
    private:
    ArmatusParams * params;

    class SubProblem {
        public:
        double score;
        size_t backPointer;
        size_t backOptimalIndex;
        bool operator<( const SubProblem & other ) const {
           return score < other.score;
        }
        bool operator==( const SubProblem & other) const {
           return (score == other.score and
                   backPointer == other.backPointer and
                   backOptimalIndex == other.backOptimalIndex);
        }
    };

    public:

    using SubProbMatrix = vector<vector<SubProblem>>;
    SubProbMatrix OPT;
    SubProbMatrix OPTD;

    explicit ArmatusDAG(ArmatusParams &p);
    
    double q(size_t k, size_t l);
    
    double s(size_t k, size_t l);

    void build();

    void computeTopK();

    DomainSet extractDomains(size_t i);
    WeightedDomainEnsemble extractTopK();
};

#endif // __ARMATUS_DAG_HPP__
