/* 
	Authors: Darya Filippova, Geet Duggal, Rob Patro
	dfilippo | geet | robp @cs.cmu.edu

	01 July 2013

	Parses 3C matrix in the format of:

		chr<chromosome_number>	<start bp>	<end bp>	<frequency1>	<frequency2> ....\n

*/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <memory>

#include <boost/range/irange.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "ArmatusUtil.hpp"
#include "ArmatusParams.hpp"
#include "IntervalScheduling.hpp"

using namespace std;

shared_ptr<SparseMatrix> parseGZipMatrix(string path) {
	cout << "parsing gzip " << path << endl;

    auto m = make_shared<SparseMatrix>();

	ifstream file(path, ios_base::in | ios_base::binary);
    assert(file.good());
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);

    string line;
    std::istream incoming(&in);
    bool firstLine = true;
    size_t i = 0;
    double tot = 0.0;
    size_t nedge = 0;
    while ( getline(incoming, line) ) {
    	vector<string> parts;
        boost::trim(line);
		boost::split(parts, line, boost::is_any_of("\t"));
        
        if (firstLine) { 
            m->resize(parts.size() - 3, parts.size() - 3, false);
            firstLine = false;
        }
        
        for (size_t j : boost::irange(3 + i, parts.size())) {
            double e  = stod(parts[j]);
            if (e > 0.0) {
                size_t row = j - 3;
                m->push_back(i, row, e);
                tot += e;
                nedge++;
            }   
        }

    	++i;
        if ( i % 1000 == 0 ) { std::cerr << "line " << i << "\n"; }
    }
    std::cerr << "Avg edge " << tot / nedge << " sum " << tot  << " cnt " << nedge << "\n";
    return m;
	// SparseSymmetricMatrix symMat(m);
 //    std::cerr << "M is " << symMat.size1() << " x " << symMat.size2() << ", with " << m.nnz() << " non-zero entries\n";
}

Domain::Domain(size_t s, size_t e) : start(s), end(e) { }

double Domain::score(ArmatusParams& p) {
    size_t d = end-start;
    return std::max((p.sums(start, end)/ std::pow(static_cast<double>(d),p.gamma)) - p.mu[d], 0.0);
}

DomainSet consensus(DomainEnsemble dEnsemble) {
    using PersistenceMap = map<Domain, int>;
    PersistenceMap pmap;

    for (auto dSet : dEnsemble) {
        for (auto domain : dSet) {
            if ( pmap.find(domain) == pmap.end() ) pmap[domain] = 0;
            pmap[domain]++;
        }
    }

    Intervals ivals;

    for (auto domainPersistence : pmap) {
        auto domain = domainPersistence.first;
        auto persistence = domainPersistence.second;
        ivals.push_back(WeightedInterval(domain.start, domain.end, persistence));      
    }

    IntervalScheduler scheduler(ivals);
    scheduler.computeSchedule();

    DomainSet dSet;
    for (auto ival : scheduler.extractIntervals()) {
        dSet.push_back(Domain(ival.start, ival.end));
    }

    return dSet;
}
