/* 
	Authors: Darya Filippova, Geet Duggal, Rob Patro
	dfilippo | geet | robp @cs.cmu.edu
        See LICENSE.txt included with this distribution.
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include <boost/range/irange.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include "ArmatusUtil.hpp"
#include "ArmatusParams.hpp"
#include "IntervalScheduling.hpp"
#include "ArmatusParams.hpp"
#include "ArmatusDAG.hpp"


MatrixProperties parseGZipMatrix(string path, int resolution, string chrom) {
	MatrixProperties prop;

    prop.matrix = std::make_shared<SparseMatrix>();

	ifstream file(path, ios_base::in | ios_base::binary);
    if (!file.good()) {
        std::cerr << "Couldn't read file " << path << std::endl;
        std::exit(1);
    }
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
            //if (parts.size() - 3 < 100) {
            if (parts.size() < 100) {
                cerr << "[INFO] Matrix is smaller than the recommended minimum size of 101." << endl;
                // exit(1);
            }
            prop.matrix->resize(parts.size(), parts.size(), false);
            prop.chrom = chrom;
            prop.resolution = resolution;
            cerr << prop.chrom << " at resolution " << prop.resolution << "bp" << endl;
            firstLine = false;
        }

        if (parts.size() != prop.matrix->size2()) {
            std::cerr << "Error: row " << i << " has " << parts.size() 
                      << " entries, but I was expecting " << prop.matrix->size2() << std::endl;
            std::exit(3);
        }

        if (i >= prop.matrix->size2()) {
            std::cerr << "Error: I was expecting " << prop.matrix->size2() 
                      << " rows, but there are more than that in the matrix file." << std::endl;
            std::exit(3);
        }

        
        for (size_t j : boost::irange(i, parts.size())) {
            double e  = stod(parts[j]);
            if (e > 0.0) {
                size_t row = j;
                prop.matrix->insert_element(i, row, e);
                tot += e;
                nedge++;
            }   
        }

    	++i;
        if ( i % 1000 == 0 ) { std::cerr << "line " << i << "\n"; }
        if (incoming.eof()) break;
    }

    // check that we read # of rows = to the number of columns in the matrix
    if (prop.matrix->size2() != i) {
        std::cerr << "Error: it doesn't look like your matrix file had enough rows." << std::endl;
        std::cerr << "Error: expecting " << prop.matrix->size2() << " but saw " << i << std::endl;
        std::exit(1);
    }
    return prop;
}

MatrixProperties parseRaoMatrix(string path, int resolution, string chrom, bool noNormalization) {
	MatrixProperties prop;
    int n;
    int M = 0;
    auto rawCountPath = path + ".RAWobserved";

    {
	    ifstream file(rawCountPath);
        if (!file.good()) {
            std::cerr << "Couldn't read file: " << rawCountPath << std::endl;
            std::exit(1);
        }
        assert(file.good());

        string line;
        int v, w;
        double count;
        int maxv = 0;
        while ( file >> v >> w >> count )  {
            if (v > maxv) maxv = v;
            if (w > maxv) maxv = w;
            M += 1;
        }
        n = maxv/resolution+1;
        prop.chrom = chrom;
        prop.resolution = resolution;
        cerr << "Building matrix for chromosome " << prop.chrom << " at resolution " << prop.resolution << "bp with " << n << " rows." << endl;
    }

    vector<double> KRnormalization;
    auto KRPath = path + ".KRnorm";
    if (!noNormalization)
    {
        cerr << "Reading KR normalization counts" << endl;
	    ifstream file(KRPath);
        string line;
        while( getline(file, line) ) {
            KRnormalization.push_back(stod(line));
        }
    }
    cout << "Initializing matrix to zero elements" << endl;
    prop.matrix = std::make_shared<SparseMatrix>(n,n,M);
    for (int i = 0; i < n; i++) {
        for (int j =0; j < n; j++) {
            prop.matrix->insert_element(i,j, 0.0);
        }
    }
    {
	    ifstream file(rawCountPath);

        string line;
        int v, w;
        double count;
        int m = 0;
        while ( file >> v >> w >> count )  {
            size_t i = v/resolution;
            size_t j = w/resolution;
            if (noNormalization) {
                prop.matrix->insert_element(i, j, log(count));
            }
            else {
                double ki = KRnormalization[i];
                double kj = KRnormalization[j];
                if (!std::isnan(ki) and !std::isnan(kj) and ki > 0 and kj > 0) {
                    prop.matrix->insert_element(i, j, log(count/(ki*kj)));
                    // cout << i << "\t" << j << "\t" << count << "\t" << count/(ki*kj) << endl;
                }
                else {
                    // Maybe write NaN or inf issue to log/warning file
                }
            }
            m++;
            if ( m % 100000 == 0 ) { std::cerr << "" << float(m)/M*100 << "%\n"; }
        }
    }

    return prop;
}

// domain size
double d(size_t const & i, size_t const & j) {return j - i + 1;}

Domain::Domain(size_t s, size_t e) : start(s), end(e) { }

DomainSet consensusDomains(WeightedDomainEnsemble& dEnsemble) {
    using PersistenceMap = map<Domain, double>;
    PersistenceMap pmap;

    for (auto dSetIdx : boost::irange(size_t{0}, dEnsemble.domainSets.size())) {
        auto& dSet = dEnsemble.domainSets[dSetIdx];
        auto weight = dEnsemble.weights[dSetIdx];
        for (auto& domain : dSet) {
            if ( pmap.find(domain) == pmap.end() ) pmap[domain] = 0;
            pmap[domain] += weight;
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

    sort(dSet.begin(), dSet.end());

    return dSet;
}

WeightedDomainEnsemble multiscaleDomains(std::shared_ptr<SparseMatrix> A, 
    float gammaMax, double stepSize, int k, int minMeanSamples, bool justThisGamma) {

    WeightedDomainEnsemble dEnsemble;
    double eps = 1e-5;
    double gamma =0.0;
    if (justThisGamma)  {
        gamma = gammaMax;
    }

    for (; gamma <= gammaMax+eps; gamma+=stepSize) {

        cerr << "gamma=" << gamma << endl;
 
        ArmatusParams params(A, gamma, k, minMeanSamples); // k parameter is not used for anything in Params
        ArmatusDAG G(params); // but is used in the DAG
        G.build();
        G.computeTopK();

        auto domainEnsemble = G.extractTopK();
        auto& domains = domainEnsemble.domainSets;
        auto& weights = domainEnsemble.weights;
        dEnsemble.domainSets.insert(dEnsemble.domainSets.end(), domains.begin(), domains.end());
        dEnsemble.weights.insert(dEnsemble.weights.end(), weights.begin(), weights.end());
        for (int i = 0; i<k; i++) {
            dEnsemble.resolutions.push_back(gamma);
            dEnsemble.optidx.push_back(i);
        }
    }

    return dEnsemble;
}

void outputDomains(DomainSet dSet, string fname, MatrixProperties matProp) {
    ofstream file;
    file.open(fname);
    int res = matProp.resolution;
    for (auto d : dSet) {
        file << matProp.chrom << "\t" << (d.start)*res << "\t" << (d.end+1)*res-1 << endl;
    }
    file.close();
}

void sanityCheck(WeightedDomainEnsemble e) {
    for (auto dset : e.domainSets) {
    }
}



