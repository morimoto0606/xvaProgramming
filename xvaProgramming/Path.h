#ifndef PATH_H_INCLUDED
#define PATH_H_INCLUDED

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include "Dual.h"
#include "Functions.h"
#include "function_traits.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;
	template <typename T>
	class Path {
	public:
		typedef T value_type;
		typedef std::size_t size_type;
		Path(const value_type& x0, const value_type& mu,
			const value_type& sigma, const size_type pathNum,
			const size_type gridNum, const double dt,
			const size_type seed)
		{
			//Generate Random Number by MT
			boost::mt19937 gen(seed);
			boost::normal_distribution<double> 
				dist(0.0, 1.0);
			boost::variate_generator <
				boost::mt19937,
				boost::normal_distribution<double >>
				rand(gen, dist);

			//Generate dB(t)
			ublas::matrix<double> dBm(pathNum, gridNum);
			for (std::size_t i = 0; i < pathNum; ++i) {
				for (std::size_t j = 0; j < gridNum; ++j) {
					dBm(i, j) = std::sqrt(dt) * rand();
				}
			}
			//Generate path of log(X(t)) by EM 
			_pathMatrix.resize(pathNum, gridNum + 1);
			for (std::size_t i = 0; i < pathNum; ++i) {
				_pathMatrix(i, 0) = cva::log(x0);
			}
			for (std::size_t i = 0; i < pathNum; ++i) {
				for (std::size_t j = 1; j < gridNum + 1; ++j) {
					_pathMatrix(i, j) = _pathMatrix(i, j - 1) 
						+ (mu - sigma * sigma / 2.0) * dt 
						+ sigma * dBm(i, j - 1);
				}
			}
			//transform from log(X(t)) to X(t)
			for (std::size_t i = 0; i < pathNum; ++i) {
				for (std::size_t j = 0; j < gridNum + 1; ++j) {
					_pathMatrix(i, j) 
						= cva::exp(_pathMatrix(i, j));
				}
			}
		}

		value_type getPathValue(
			const size_type pathIndex,
			const size_type gridIndex) const
		{
			return _pathMatrix(pathIndex, gridIndex);
		}

		ublas::vector<value_type>
			getTimewisePath(const size_type i) const
		{
			return ublas::row(_pathMatrix, i);
		}
		ublas::vector<value_type>
			getPathwisePath(const size_type j) const
		{
			return ublas::column(_pathMatrix, j);
		}
		const ublas::matrix<value_type>&
			getPathMatrix() const
		{
			return _pathMatrix;
		}
		size_type pathNum() const {
			return _pathMatrix.size1();
		}
		size_type gridNum() const {
			return _pathMatrix.size2() - 1;
		}
	private:
		ublas::matrix<value_type> _pathMatrix;
	};
} //namespace cva
#endif