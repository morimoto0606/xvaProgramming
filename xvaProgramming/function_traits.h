#ifndef FUNCTION_TRAITS_H_INCLUDED
#define FUNCTION_TRAITS_H_INCLUDED
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include "Dual.h"

namespace cva {
	template <typename T>
	struct zero_floor_traits {
		typedef T value_type;
		typedef T result_type;
		static const result_type
			apply(const value_type& x)
		{
			return std::max(x, 0.0);
		}
	};

	template <typename T>
	struct zero_floor_traits<Dual<T> > {
		typedef Dual<T> value_type;
		typedef Dual<T> result_type;
		static const result_type
			apply(const value_type& x)
		{
			const T value = std::max(x.value(), 0.0);
			const T deriv = x.value() >= 0 ? x.deriv() : 0.0;
			return result_type(value, deriv);
		}
	};

	template <typename T>
	struct exp_traits {
	public:
		typedef T value_type;
		typedef T result_type;
		static const result_type apply(const value_type&x)
		{
			return std::exp(x);
		}
	};

	template <typename T>
	struct exp_traits<Dual<T> > {
	public:
		typedef Dual<T> value_type;
		typedef Dual<T> result_type;
		static const result_type apply(const value_type& x)
		{
			const T value = std::exp(x.value());
			const T deriv = x.deriv() * value;
			return result_type(value, deriv);
		}
	};

	template <typename T>
	struct log_traits {
		typedef T value_type;
		typedef T result_type;
		static const result_type apply(const value_type& x)
		{
			return std::log(x);
		}
	};

	template <typename T>
	struct log_traits<Dual<T> > {
		typedef Dual<T> value_type;
		typedef Dual<T> result_type;
		static const result_type apply(const value_type& x)
		{
			return result_type(std::log(x.value()), x.deriv() / x.value());
		}
	};

	template <typename T>
	struct normal_cdf_traits {
		typedef T value_type;
		typedef T result_type;
		static const result_type apply(const value_type& x)
		{
			boost::math::normal_distribution<> normal;
			return boost::math::cdf(normal, x);
		}
	};

	template <typename T>
	struct normal_cdf_traits<Dual<T> > {
		typedef Dual<T> value_type;
		typedef Dual<T> result_type;
		static const result_type apply(const value_type& x)
		{
			boost::math::normal_distribution<> normal;
			return result_type(boost::math::cdf(normal, x.value()),
				x.deriv() * boost::math::pdf(normal, x.value()));
		}
	};
	
	template <typename T>
	struct sqrt_traits {
		typedef T value_type;
		typedef T result_type;
		static const result_type apply(const value_type& x)
		{
			return std::sqrt(x);
		}
	};
	template <typename T>
	struct sqrt_traits<Dual<T> > {
		typedef Dual<T> value_type;
		typedef Dual<T> result_type;
		static const result_type apply(const value_type& x)
		{
			const T value = std::sqrt(x.value());
			const T deriv = 0.5 / value * x.deriv();
			return Dual<double>(value, deriv);
		}
	};


} //namespace cva
#endif