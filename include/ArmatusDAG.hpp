/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
*/

#ifndef __ARMATUS_DAG_HPP__
#define __ARMATUS_DAG_HPP__

#include "ArmatusUtil.hpp"
#include "ArmatusParams.hpp"

using namespace boost;
using namespace std;

class ArmatusDAG {
    private:
    ArmatusParams * params;

    class BackPointer {
        public:
        size_t edge;
        size_t childSolution;
        double score;
        bool operator==(const BackPointer &other) const { 
            return (edge == other.edge) and (childSolution == other.childSolution) and
                   (score == other.score);
        }
        bool operator!=(const BackPointer &other) const { 
            return !((*this) == other);
        }
    };

    class SubProblem {
        public:
        vector<BackPointer> topK;
    };

    using SubProbsVec = vector<SubProblem>;
    SubProbsVec subProbs;

    using EdgeWeightsVec = vector< vector<double> >;
    EdgeWeightsVec edgeWeights;

    public:
    
    explicit ArmatusDAG(ArmatusParams &p);
    
    void build();

    double q(size_t k, size_t l);
    
    double s(size_t k, size_t l);

    DomainSet viterbiPath();

    void computeTopK(uint32_t k);
};

#endif // __ARMATUS_DAG_HPP__
