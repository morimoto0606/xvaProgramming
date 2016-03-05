#ifndef PATHMAKER_H_INCLUDED
#define PATHMAKER_H_INCLUDED

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
#include "Path.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;
	template <typename T>
	struct pathmaker_traits {
	public:
		typedef T value_type;
		typedef Path<value_type> result_type;
		static const value_type apply(
			const double value, 
			const double deriv)
		{
			return value_type(value, deriv);
		}
	};
	template <>
	struct pathmaker_traits<double> {
		typedef double value_type;
		typedef Path<double> result_type;
		static const value_type apply(
			const double value, const double deriv)
		{
			return value;
		}
	};

	enum shockTypeEnum { undEnum = 0, volEnum = 1 };
	enum productTypeEnum {
		fwdEnum = 0, eurEnum = 1,
		rrEnum = 2, mountainEnum = 3
	};

	template <typename T>
	typename pathmaker_traits<T>::result_type makePath(
		const double x0, const double mu0, const double sigma0,
		const double dt, const std::size_t gridNum,
		const std::size_t pathNum,
		const shockTypeEnum shockType, const std::size_t seed)
	{
		const T x = shockType == undEnum
			? pathmaker_traits<T>::apply(x0, 1.0)
			: pathmaker_traits<T>::apply(x0, 0.0);
		const T sigma = shockType == volEnum
			? pathmaker_traits<T>::apply(sigma0, 1.0)
			: pathmaker_traits<T>::apply(sigma0, 0.0);
		const T mu(mu0);
		return pathmaker_traits<T>::result_type(x, mu,
			sigma, pathNum, gridNum, dt, seed);
	}
} //namespace cva
#endif