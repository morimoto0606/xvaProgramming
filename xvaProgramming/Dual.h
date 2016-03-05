#pragma once
#include <boost/numeric/ublas/functional.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
namespace cva {
	template <typename T>
	class Dual {
	public:
		typedef T value_type;
		typedef T result_type;
		//constructors
		Dual() : _value(0.0), _deriv(0.0) {}

		explicit Dual(const value_type& value)
			: _value(value), _deriv(0.0) {}

		Dual(const value_type& value, const value_type& deriv)
			: _value(value), _deriv(deriv) {}

		Dual(const Dual& rhs) {
			_value = rhs._value;
			_deriv = rhs._deriv;
		}

		Dual& operator =(const Dual<value_type>& rhs) {
			if (this != &rhs) {
				_value = rhs._value;
				_deriv = rhs._deriv;
			}
			return *this;
		}

		//accessors
		const result_type& value() const {
			return _value;
		}
		const result_type& deriv() const {
			return _deriv;
		}

		result_type operator()(const Dual<value_type>& rhs) const
		{
			return rhs.value();
		}
		Dual& operator+=(const Dual<value_type>& rhs)
		{
			_value += rhs._value;
			_deriv += rhs._deriv;
			return *this;
		}
		Dual& operator+=(const value_type& rhs)
		{
			_value += rhs;
			return *this;
		}

		Dual& operator-=(const Dual<value_type>& rhs)
		{
			_value -= rhs._value;
			_deriv -= rhs._deriv;
			return *this;
		}
		Dual& operator-=(const value_type& rhs)
		{
			_value -= rhs;
			return *this;
		}
		Dual& operator*=(const Dual<value_type>& rhs)
		{
			const value_type value1 = _value;
			const value_type deriv1 = _deriv;
			const value_type value2 = rhs._value;
			const value_type deriv2 = rhs._deriv;
			_value = value1 * value2;
			_deriv = deriv1 * value2 + value1 * deriv2;
			return *this;
		}
		Dual& operator*=(const value_type& rhs)
		{
			_value *= rhs;
			_deriv *= rhs;
			return *this;
		}

		Dual& operator/=(const Dual<value_type>& rhs)
		{
			const value_type value1 = _value;
			const value_type deriv1 = _deriv;
			const value_type value2 = rhs._value;
			const value_type deriv2 = rhs._deriv;
			_value = value1 / value2;
			_deriv = (deriv1 * value2 - value1 * deriv2) / (value2 * value2);
			return *this;
		}
		Dual& operator/=(const value_type& rhs)
		{
			_value /= rhs;
			_deriv /= rhs;
			return *this;
		}
	private:
		value_type _value;
		value_type _deriv;
	};

	//operators
	template <typename T>
	Dual<T> operator +(const Dual<T>& lhs, const Dual<T>& rhs)
	{
		return Dual<T>(lhs.value() + rhs.value(),
			lhs.deriv() + rhs.deriv());
	}
	template <typename T>
	Dual<T> operator +(const Dual<T>& lhs, const T& rhs)
	{
		return Dual<T>(lhs.value() + rhs, lhs.deriv());
	}
	template <typename T>
	Dual<T> operator +(const T& lhs, const Dual<T>& rhs)
	{
		return rhs + lhs;
	}
	template <typename T>
	Dual<T> operator -(const Dual<T>& lhs, const Dual<T>& rhs)
	{
		return Dual<T>(lhs.value() - rhs.value(),
			lhs.deriv() - rhs.deriv());
	}
	template <typename T>
	Dual<T> operator -(const Dual<T>& lhs, const T& rhs)
	{
		return Dual<T>(lhs.value() - rhs, lhs.deriv());
	}
	template <typename T>
	Dual<T> operator -(const T& lhs, const Dual<T>& rhs)
	{
		return Dual<T>(lhs - rhs.value(), -lhs.deriv());
	}

	template <typename T>
	Dual<T> operator *(const Dual<T>& lhs, const Dual<T>& rhs)
	{
		const T lvalue = lhs.value();
		const T rvalue = rhs.value();
		const T lderiv = lhs.deriv();
		const T rderiv = rhs.deriv();

		return Dual<T>(lvalue * rvalue, lderiv * rvalue + lvalue * rderiv);
	}
	template <typename T>
	Dual<T> operator *(const Dual<T>& lhs, const T& rhs)
	{
		return Dual<T>(lhs.value() * rhs, lhs.deriv() * rhs);
	}
	template <typename T>
	Dual<T> operator *(const T& lhs, const Dual<T>& rhs)
	{
		return rhs * lhs;
	}
	template <typename T>
	Dual<T> operator /(const Dual<T>&lhs, const Dual<T>& rhs)
	{
		const T lvalue = lhs.value();
		const T rvalue = rhs.value();
		const T lderiv = lhs.deriv();
		const T rderiv = rhs.deriv();
		return Dual<T>(lvalue / rvalue,
			(lderiv * rvalue - lvalue * rderiv) / (rvalue * rvalue));
	}
	template <typename T>
	Dual<T> operator /(const Dual<T>&lhs, const T& rhs)
	{
		return Dual<T>(lhs.value() / rhs, lhs.deriv() / rhs);
	}
	template <typename T>
	Dual<T> operator /(const T&lhs, const Dual<T>& rhs)
	{
		return Dual<T>(lhs.value() / rhs,
			-lhs * rhs.deriv() / (rhs.value() * rhs.value()));
	}
	template <typename T>
	bool operator == (const Dual<T>& lhs, const Dual<T>& rhs) {
		return lhs.value() == rhs.value() && lhs.deriv() == rhs.deriv();
	}
	template <typename T>
	bool operator == (const Dual<T>& lhs, const T& rhs) {
		return lhs.value() == rhs;
	}
	template <typename T>
	bool operator == (T& lhs, const Dual<T>& rhs) {
		return rhs == lhs;
	}
	template <typename T>
	bool operator != (const Dual<T>& lhs, const Dual<T>& rhs) {
		return !(lhs == rhs);
	}
	template <typename T>
	bool operator != (const Dual<T>& lhs, const T& rhs) {
		return !(lhs == rhs);
	}
	template <typename T>
	bool operator != (const T& lhs, const Dual<T>& rhs) {
		return !(lhs == rhs);
	}
	template <typename T>
	bool operator < (const Dual<T>& lhs, const Dual<T>& rhs) {
		return lhs.value() < rhs.value();
	}
	template <typename T>
	bool operator < (const Dual<T>& lhs, const T& rhs) {
		return lhs.value() < rhs;
	}
	template <typename T>
	bool operator < (const T& lhs, const Dual<T>& rhs) {
		return lhs < rhs.value();
	}
	template <typename T>
	bool operator > (const Dual<T>& lhs, const Dual<T>& rhs) {
		return lhs.value() > rhs.value();
	}
	template <typename T>
	bool operator > (const Dual<T>& lhs, const T& rhs) {
		return lhs.value() > rhs;
	}
	template <typename T>
	bool operator > (const T& lhs, const Dual<T>& rhs) {
		return lhs > rhs.value();
	}

} //namespace cva

namespace boost {
	namespace numeric {
		namespace ublas {
			template<class F>
			struct matrix_scalar_unary_traits<cva::Dual<double>, F> {
				typedef matrix_scalar_unary<double, F> expression_type;
#ifndef BOOST_UBLAS_SIMPLE_ET_DEBUG
				typedef double result_type;
#else
				typedef typename double result_type;
#endif
			};
		}
	}
}
