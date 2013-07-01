/* 
	Authors: Darya Filippova, Geet Duggal, Rob Patro
	dfilippo | geet | robp @cs.cmu.edu

	01 July 2013

	Parses 3C matrix in the format of:

		chr<chromosome_number>	<start bp>	<end bp>	<frequency1>	<frequency2> ....\n

*/

#include "MatrixParser.hpp"
#include <iostream>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;

void MatrixParser::parseGZipMatrix(string path) {
	cout << "parsing gzip " << path << endl;
}