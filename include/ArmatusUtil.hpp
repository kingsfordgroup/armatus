#ifndef __ARMATUS_UTIL_HPP__
#define __ARMATUS_UTIL_HPP__

#include <memory>
#include <string>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using SparseMatrix = boost::numeric::ublas::compressed_matrix<double>;

std::shared_ptr<SparseMatrix> parseGZipMatrix(std::string path);

#endif // __ARMATUS_UTIL_HPP__