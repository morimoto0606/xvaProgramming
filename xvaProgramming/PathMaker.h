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
#include "promote_traits.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;
	//template <typename T>
	//struct pathmaker_traits {
	//public:
	//	typedef T value_type;
	//public:
	//	static const Path<value_type> apply(
	//		const double value, 
	//		const double deriv)
	//	{
	//		return value_type(value, deriv);
	//	}
	//};
	//template <>
	//struct pathmaker_traits<double> {
	//	typedef double value_type;
	//	typedef Path<double> result_type;
	//	static const value_type apply(
	//		const double value, const double deriv)
	//	{
	//		return value;
	//	}
	//};

	enum shockTypeEnum { undEnum = 0, volEnum = 1 };
	enum productTypeEnum {
		fwdEnum = 0, eurEnum = 1,
		rrEnum = 2, mountainEnum = 3
	};

	//template <typename X,  typename S>
	//Path<typename promote_traits<X, S>::result_type> makePath(
	//	const X& x0, const double mu0, const S& sigma0,
	//	const double dt, const std::size_t gridNum,
	//	const std::size_t pathNum,
	//	const std::size_t seed)
	//{
	//	return Path<typename promote_traits<X, S>::result_type>(x0, mu0,
	//		sigma0, pathNum, gridNum, dt, seed);
	//}
} //namespace cva
#endif