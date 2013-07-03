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
#include "MatrixParser.hpp"

using namespace std;

shared_ptr<MatrixParser::SparseMatrix> MatrixParser::parseGZipMatrix(string path) {
	cout << "parsing gzip " << path << endl;

    auto m = make_shared<SparseMatrix>();

	ifstream file(path, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
	boost::numeric::ublas::mapped_matrix<double> m (3, 3, 3 * 3);

    string line;
    std::istream incoming(&in);
    bool firstLine = true;
    size_t i = 0;
    while ( getline(incoming, line) ) {
    	vector<string> parts;
		boost::split(parts, line, boost::is_any_of("\t"));
        
        if (firstLine) { 
            m->resize(parts.size() - 3, parts.size() - 3, false);
            firstLine = false;
        }
        
        for (size_t j : boost::irange(3 + i, parts.size())) {
            float e  = stof(parts[j]);
            if (e > 0.0) {
                size_t row = j - 3;
                m->push_back(i, row, e);
            }
        }

    	++i;
        if ( i % 1000 == 0 ) { std::cerr << "line " << i << "\n"; }
    }

    return m;
	// SparseSymmetricMatrix symMat(m);
 //    std::cerr << "M is " << symMat.size1() << " x " << symMat.size2() << ", with " << m.nnz() << " non-zero entries\n";
}