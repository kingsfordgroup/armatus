#ifndef __ARMATUS_DAG_HPP__
#define __ARMATUS_DAG_HPP__


#include "ArmatusParams.hpp"

using namespace boost;
using namespace std;

class ArmatusDAG {
    private:
    ArmatusParams * params;

    class BackPointer {
        public:
        EdgeID edge;
        size_t childSolution;
        double score;
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
    
    ArmatusDAG(ArmatusParams &p) : 
        solutionDAG(SubProbsVec(p.n+1)), 
        edgeWeights(EdgeWeightsVec(p.n+1, vector<double>())) {
        params = &p;
    }
    
    void build() {
        for (size_t i = 1; i <= params->n; i++) {
            edgeWeights[i]
        }
    }
};

#endif // __ARMATUS_DAG_HPP__
