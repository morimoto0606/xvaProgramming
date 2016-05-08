#ifndef FUNCTION_TRAITS_H_INCLUDED
#define FUNCTION_TRAITS_H_INCLUDED
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include "Dual.h"
#include <algorithm>
#include <vector>
#include <iterator>

namespace cva {
	namespace ublas = boost::numeric::ublas;

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
			return std::log(std::abs(x));
		}
	};

	template <typename T>
	struct log_traits<Dual<T> > {
		typedef Dual<T> value_type;
		typedef Dual<T> result_type;
		static const result_type apply(const value_type& x)
		{
			return result_type(std::log(std::abs(x.value())), x.deriv() / x.value());
		}
	};

	template <typename T>
	struct power_traits {
		typedef T value_type;
		typedef T result_type;
		static const result_type
			apply(const value_type& x, const std::size_t d)
		{
			const T sign = (x < 0.0) ? std::pow(-1.0, d) : 1.0;
			return sign * std::pow(std::abs(x), d);
		}
	};

	template <typename T>
	struct  power_traits<Dual<T>>
	{
		typedef Dual<T> value_type;
		typedef Dual<T> result_type;
		static const result_type
			apply(const value_type& x, const std::size_t d)
		{
			const T sign = (x.value() < 0.0) ? std::pow(-1.0, d) : 1.0;
			const value_type y(std::abs(x.value()), std::abs(x.deriv()));
			return sign * exp_traits<Dual<double>>::apply(
				static_cast<double>(d) * log_traits<Dual<double>>::apply(y));
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


	template <typename T>
	struct average_traits {
		template<typename It>
		static const typename std::iterator_traits<It>::value_type
			apply(const It itFirst, const It itLast,
				const typename std::iterator_traits<It>::value_type& ini)
		{
			return std::accumulate(itFirst, itLast, ini) / static_cast<double>(itLast - itFirst);
		}
	};

	template <typename T>
	struct average_traits <Dual<T>>{
		template<typename It>
		static const typename std::iterator_traits<It>::value_type
			apply(const It itFirst, const It itLast, 
				const typename std::iterator_traits<It>::value_type& ini)
		{
			std::vector<T> value;
			std::vector<T> deriv;
			std::for_each(itFirst, itLast, [&value](const Dual<T>& dual) mutable->void {value.push_back(dual.value()); });
			std::for_each(itFirst, itLast, [&deriv](const Dual<T>& dual) mutable->void {deriv.push_back(dual.deriv()); });
			T v = std::accumulate(value.begin(), value.end(), 0.0);
			T d = std::accumulate(deriv.begin(), deriv.end(), 0.0);
			return typename std::iterator_traits<It>::value_type(v, d) / static_cast<double>(itLast - itFirst);
		}
	};

} //namespace cva
#endif