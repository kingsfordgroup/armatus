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
        vector<double> score;
    };

    vector<SubProblem> OPT;
    vector<SubProblem> OPTD;

    public:
    
    explicit ArmatusDAG(ArmatusParams &p);
    
    double q(size_t k, size_t l);
    
    double s(size_t k, size_t l);

    void build();

    void computeTopK(uint32_t k);

    WeightedDomainEnsemble extractTopK(uint32_t k);
};

#endif // __ARMATUS_DAG_HPP__
