/* 
	Authors: Darya Filippova, Geet Duggal, Rob Patro
	dfilippo | geet | robp @cs.cmu.edu

	01 July 2013

	Parses 3C matrix in the format of:

		chr<chromosome_number>	<start bp>	<end bp>	<frequency1>	<frequency2> ....\n

*/

#include "MatrixParser.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

void MatrixParser::parseGZipMatrix(string path) {
	cout << "parsing gzip " << path << endl;

	ifstream file(path, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    // boost::iostreams::copy(in, cout);

    string line;
    std::istream incoming(&in);
    int i = 0;
    while ( incoming ) {
    	if (i > 10) break;
    	getline(incoming, line);
    	cout << line << endl;
    	vector<string> parts;
		// parts.reserve(3);
		boost::split(parts, line, boost::is_any_of("\t"));
		cout << parts.size() << endl;
    	i++;
    }
	
}