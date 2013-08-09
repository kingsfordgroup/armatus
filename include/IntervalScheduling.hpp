/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
*/

#ifndef  __INTERVAL_SCHEDULING_HPP__
#define  __INTERVAL_SCHEDULING_HPP__

#include <vector>

using namespace std;

class WeightedInterval {
    public:
    int start;
    int end;
    double score;

    WeightedInterval(int s, int e, double sc);

    bool operator< (const WeightedInterval &other) const {
        return end < other.end;
    }
};

using Intervals = vector<WeightedInterval>;

class IntervalScheduler {
    private:
    size_t n;
    Intervals ivals;
    vector<double> bestGeneScore;
    vector<size_t> p;
    size_t previousDisjointInterval(size_t j);

    public:
    IntervalScheduler(Intervals &inputIntervals);
    void computeSchedule();
    Intervals extractIntervals();
};

#endif // __INTERVAL_SCHEDULING_HPP__