/* 
    Authors: Darya Filippova, Geet Duggal, Rob Patro
    dfilippo | geet | robp @cs.cmu.edu
    See LICENSE.txt included with this distribution.
*/


#include "IntervalScheduling.hpp"
#include <algorithm>

WeightedInterval::WeightedInterval(int s, int e, double sc) {
    if (s < e) {
        start = s;
        end = e;
    } else {
        start = e;
        end = s;
    }
    score = sc;
}

size_t IntervalScheduler::previousDisjointInterval(size_t j) {
    for (size_t i = j; i > 0; i--) {
        if (ivals[i].end < ivals[j].start) return i;
    }
    return 0;
}


IntervalScheduler::IntervalScheduler(Intervals &inputIntervals) : 
    n(inputIntervals.size()),
    bestGeneScore(vector<double>(n+1, 0)),
    p(vector<size_t>(n+1, 0)) { 
        
        sort(inputIntervals.begin(), inputIntervals.end());

        // We want the gene list to start at index 1
        ivals.push_back(WeightedInterval(-1,-1,-1));    
        for (auto i : inputIntervals) ivals.push_back(i);

        for (size_t j = 1; j < n+1; j++) {
            p[j] = previousDisjointInterval(j);
        }
}

void IntervalScheduler::computeSchedule() {

    for (size_t j = 1; j < n+1; j++) {
        double scoreIfChosen = ivals[j].score + bestGeneScore[p[j]];
        double scoreIfIgnored = bestGeneScore[j-1];

        bestGeneScore[j] = max(scoreIfChosen, scoreIfIgnored);
    } 
}

Intervals IntervalScheduler::extractIntervals() {
    Intervals extracted;
    size_t j = bestGeneScore.size()-1;
    while (j > 0) {
        if (bestGeneScore[j] != bestGeneScore[j-1]) {
            extracted.push_back(ivals[j]);
            j = previousDisjointInterval(j);
        } else { j--; }
    }
    return extracted;
}
