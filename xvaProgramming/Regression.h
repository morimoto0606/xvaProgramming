#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/function.hpp>
#include "PayOff.h"
#include "Path.h"
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp>     
#include <boost/numeric/ublas/io.hpp>
#include <boost/type_traits.hpp>
#include "BasisFunctions.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;

	template <typename T>
	ublas::matrix<T> getBasisMatrix(
		const T& pathValue,
		const ublas::vector<boost::function<T(const T&)> >& functions)
	{
		const std::size_t basisNum = functions.size();
		ublas::vector<T> coeffient(basisNum);
		ublas::matrix<T> basisMatrix(basisNum, basisNum);
		for (std::size_t i = 0; i < basisNum; ++i) {
			for (std::size_t j = 0; j < basisNum; ++j) {
				basisMatrix(i, j)
					= (functions(i))(pathValue) * (functions(j))(pathValue);
			}
		}
		return basisMatrix;
	}

	template <typename T, typename D>
	ublas::matrix<typename BasisFunctions<D>::state_type> getBasisMatrix(
		const std::size_t pathIndex,
		const std::size_t gridIndex,
		const Path<T>& path,
		const BasisFunctions<D>& functions)
	{
		typedef BasisFunctions<D>::state_type state_type;
		const std::size_t basisNum = functions.size();
		ublas::vector<state_type> coeffient(basisNum);
		ublas::matrix<state_type> basisMatrix(basisNum, basisNum);
		for (std::size_t i = 0; i < basisNum; ++i) {
			for (std::size_t j = 0; j < basisNum; ++j) {
				basisMatrix(i, j)
					= (functions(path, pathIndex, gridIndex))(i)
					* (functions(path, pathIndex,gridIndex))(j);
			}
		}
		return basisMatrix;
	}

	template <typename T, typename U>
	ublas::vector<T> calcPayoffMultBasis(
		const T& pathValue,
		const ublas::vector<T>& timewisePath,
		const PayOff<U>& payoff,
		const ublas::vector<boost::function<T(const T&)> >& functions)
	{
		const std::size_t basisNum = functions.size();
		T payoffValue = payoff()(timewisePath);
		ublas::vector<T> payoffMultBasis(basisNum);
		for (std::size_t i = 0; i < basisNum; ++i) {
			payoffMultBasis(i)
				= (functions(i))(pathValue) * payoffValue;
		}
		return payoffMultBasis;
	}
	
	template <typename T, typename U, typename D>
	ublas::vector<typename BasisFunctions<D>::state_type> calcPayoffMultBasis(
		const std::size_t pathIndex,
		const std::size_t gridIndex,
		const Path<T>& path,
		const PayOff<U>& payoff,
		const BasisFunctions<D>& functions)
	{
		typedef BasisFunctions<D>::state_type state_type;
		const std::size_t basisNum = functions.size();
		T payoffValue = payoff()(path.getTimewisePath(pathIndex));
		ublas::vector<state_type> payoffMultBasis(basisNum);
		for (std::size_t i = 0; i < basisNum; ++i) {
			payoffMultBasis(i)
				= functions(path, pathIndex, gridIndex)(i) * payoffValue;
		}
		return payoffMultBasis;
	}

	/* Matrix inversion routine.
	Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
	template<class T>
	bool invertMatrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse)
	{
		typedef ublas::permutation_matrix<std::size_t> pmatrix;

		// create a working copy of the input
		ublas::matrix<T> A(input);

		// create a permutation matrix for the LU-factorization
		pmatrix pm(A.size1());

		// perform LU-factorization
		int res = ublas::lu_factorize(A, pm);
		if (res != 0) {
			return false;
		}

		// create identity matrix of "inverse"
		inverse.assign(ublas::identity_matrix<T>(A.size1()));

		// backsubstitute to get the inverse
		ublas::lu_substitute(A, pm, inverse);
		return true;
	}


	template <typename T, typename U>
	ublas::vector<T> regresssion(
		std::size_t gridIndex,
		const PayOff<U>& payoff,
		const Path<T>& path,
		const ublas::vector<boost::function<T (const T&)> >& functions)
	{
		const std::size_t basisNum = functions.size();
		const std::size_t pathNum = path.pathNum();
		
		// matrix (basis(i) * basis(j))
		ublas::matrix<T> basisMatrix 
			= ublas::zero_matrix<T>(basisNum, basisNum);
		//vector (basis(i) *Payoff)
		ublas::vector<T> payoffMultBasis
			=ublas::zero_vector<T>(basisNum);

		// take expectation of basisMatrix and payoffBasis
		for (std::size_t k = 0; k< pathNum; ++k) {
			basisMatrix += getBasisMatrix(
				path.getPathValue(k, gridIndex),
				functions);
			payoffMultBasis += calcPayoffMultBasis(
				path.getPathValue(k, gridIndex),
				path.getTimewisePath(k),
				payoff, functions);
		}
		basisMatrix /= pathNum;
		payoffMultBasis /= pathNum;
		// LU Decomposition
		ublas::matrix<T> basisInverse
			=ublas::identity_matrix<T>(basisNum, basisNum);
		bool isScess = invertMatrix(basisMatrix, basisInverse);
		return ublas::prod(payoffMultBasis, basisInverse);
	}

	template <typename T, typename U, typename D>
	ublas::vector<T> regresssion(
		std::size_t gridIndex,
		const PayOff<U>& payoff,
		const Path<T>& path,
		const BasisFunctions<D>& functions)
	{
		typedef BasisFunctions<D>::state_type state_type;
		const std::size_t basisNum = functions.size();
		const std::size_t pathNum = path.pathNum();

		// matrix (basis(i) * basis(j))
		ublas::matrix<T> basisMatrix
			= ublas::zero_matrix<T>(basisNum, basisNum);
		//vector (basis(i) *Payoff)
		ublas::vector<T> payoffMultBasis
			= ublas::zero_vector<T>(basisNum);

		// take expectation of basisMatrix and payoffBasis
		for (std::size_t k = 0; k < pathNum; ++k) {
			basisMatrix += getBasisMatrix(
				k, gridIndex, path, functions);
			payoffMultBasis += calcPayoffMultBasis(
				k, gridIndex, path, payoff, functions);
		}
		basisMatrix /= pathNum;
		payoffMultBasis /= pathNum;
		// LU Decomposition
		ublas::matrix<T> basisInverse
			= ublas::identity_matrix<T>(basisNum, basisNum);
		bool isScess = invertMatrix(basisMatrix, basisInverse);
		return ublas::prod(payoffMultBasis, basisInverse);
	}
}//namespace cva {