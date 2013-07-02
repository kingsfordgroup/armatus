#ifndef __MATRIX_PARSER_HPP__
#define __MATRIX_PARSER_HPP__

#include <memory>
#include <string>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

class MatrixParser {
  using SparseMatrix = boost::numeric::ublas::compressed_matrix<double>;
  //using SparseSymmetricMatrix = ublas::sparse_adaptor<SparseMatrix, ublas::upper>;
	public: 
		std::shared_ptr<SparseMatrix> parseGZipMatrix(std::string path);
};

#endif // __MATRIX_PARSER_HPP__