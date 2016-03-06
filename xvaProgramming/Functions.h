#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "function_traits.h"
#include <boost/numeric/ublas/vector_expression.hpp>
#include <numeric>

namespace ublas = boost::numeric::ublas;
namespace cva{
	template <typename T>
	typename zero_floor_traits<T>::value_type
		zeroFloor(const T& x)
	{
		return zero_floor_traits<T>::apply(x);
	}
	
	template <typename T>
	typename exp_traits<T>::value_type
		exp(const T& x)
	{
		return exp_traits<T>::apply(x);
	}

	template <typename T>
	typename log_traits<T>::value_type
		log(const T& x)
	{
		return log_traits<T>::apply(x);
	}

	template <typename T>
	typename normal_cdf_traits<T>::value_type
		normalCdf(const T& x)
	{
		return normal_cdf_traits<T>::apply(x);
	}

	template <typename T>
	typename sqrt_traits<T>::result_type sqrt(
		const T& x)
	{
		return sqrt_traits<T>::apply(x);
	}
	//template <typename T, std::size_t n>
	//typename sqrt_traits<T>::result_type sqrt(
	//	const T& x)
	//{
	//	return sqrt_traits<T>::apply(x);
	//}

	class Monomial {
	public:
		explicit Monomial(const double order) : _order(order) {}
		
		template <typename T>
		T operator()(const T& x) const
		{
			return _order == 0.0 ? 1.0 : std::pow(x, _order);
		}
		template <typename T>
		Dual<T> operator()(const Dual<T>& x) const
		{
			T value = _order == 0.0 ? 1.0 : std::pow(x.value(), _order);
			T deriv = _order * std::pow(x.value(), _order - 1) * x.deriv();
			return Dual<T>(value, deriv);
		}
	private: 
		double _order;
	};

	class Monomial2 {
	public:
		Monomial2(const double order1, const double order2)
		: _order1(order1), _order2(order2) {}

		template <typename T>
		T operator()(const T& x, const T& y) const
		{
			return std::pow(x, _order1) * std::pow(y, _order2);
		}

		template <typename T>
		Dual<T> operator()(const Dual<T>& x, const Dual<T>& y) const
		{
			T value = std::pow(x.value(), _order1) * std::pow(y.value(), _order2);
			T deriv = _order1 * std::pow(x.value(), _order1 - 1) * x.deriv * y.value()
				+ _order2 * std::pow(y.value(), _order2 - 1) * y.deriv() * x.value();
			return Dual<T>(value, deriv);
		}
	private:
		double _order1;
		double _order2;
	};

	class PathwiseMonomial {
	public:
		PathwiseMonomial(const std::size_t grid,  std::size_t order) 
		: _grid(grid), _order(order) {}
		
		template <typename T>
		typename T::value_type
			operator()(const ublas::vector_expression<T>& x) const
		{
			return _order == 0.0 ? 1.0 : std::pow((x()(_grid)), _order);
		}
	private:
		std::size_t _order;
		std::size_t _grid;
	};

	class PathwiseSum {
	public:
		PathwiseSum(const std::size_t grid, std::size_t order)
			: _grid(grid), _order(order) {}

		template <typename T>
		typename T::value_type operator()(const ublas::vector_expression<T>& x) const
		{			
			return _order == 0.0 
				? 1.0
				: std::pow(std::accumulate(x().begin(), x().begin() + _grid + 1, 0.0), _order);
		}
	private:
		std::size_t _order;
		std::size_t _grid;
	};

} //namespace cva

namespace boost {	namespace numeric {	namespace ublas {
	// Define properties for a generic scalar type
	template<>
	struct scalar_traits<cva::Dual<double> >{
		typedef scalar_traits<cva::Dual<double> > self_type;
		typedef cva::Dual<double> value_type;
		typedef const cva::Dual<double> &const_reference;
		typedef cva::Dual<double> &reference;

		typedef double real_type;
		typedef real_type precision_type;       // we do not know what type has more precision then the real_type

		static const unsigned plus_complexity = 1;
		static const unsigned multiplies_complexity = 1;

		static
			BOOST_UBLAS_INLINE
			real_type real(const_reference t) {
			return t.value();
		}
		static
			BOOST_UBLAS_INLINE
			real_type imag(const_reference /*t*/) {
			return 0;
		}
		static
			BOOST_UBLAS_INLINE
			value_type conj(const_reference t) {
			return t;
		}

		static
			BOOST_UBLAS_INLINE
			real_type type_abs(const_reference t) {
			return sqrt(t.value() * t.value() + t.deriv() * t.deriv());
		}
		static
			BOOST_UBLAS_INLINE
			value_type type_sqrt(const_reference t) {
			// force a type conversion back to value_type for intgral types
			return value_type(cva::sqrt(t));
		}

		static
			BOOST_UBLAS_INLINE
			real_type norm_1(const_reference t) {
			return self_type::type_abs(t);
		}
		static
			BOOST_UBLAS_INLINE
			real_type norm_2(const_reference t) {
			return self_type::type_abs(t);
		}
		static
			BOOST_UBLAS_INLINE
			real_type norm_inf(const_reference t) {
			return self_type::type_abs(t);
		}

		static
			BOOST_UBLAS_INLINE
			bool equals(const_reference t1, const_reference t2) {
			return self_type::norm_inf(t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
				(std::max) ((std::max) (self_type::norm_inf(t1),
					self_type::norm_inf(t2)),
					BOOST_UBLAS_TYPE_CHECK_MIN);
		}
	};

	template<>
	struct type_traits<cva::Dual<double> >
		: scalar_traits <cva::Dual<double> > {
		typedef type_traits<cva::Dual<double> > self_type;
		typedef cva::Dual<double> value_type;
		typedef const cva::Dual<double> &const_reference;
		typedef cva::Dual<double> &reference;

		typedef double real_type;
		typedef real_type precision_type;
		static const unsigned multiplies_complexity = 1;

	};
}}}
#endif