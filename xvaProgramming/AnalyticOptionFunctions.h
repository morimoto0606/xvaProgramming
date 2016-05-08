#pragma once
#include "Functions.h"
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;
namespace cva {
	//E[aX(t) + b]
	class ForwardCalculator {
	public:
		ForwardCalculator(const double mu,
			const double gearing, const double strike, 
			const double maturity) : _mu(mu), _gearing(gearing),
			_strike(strike), _maturity(maturity) {}
		template <typename T>
		T operator()(const T& x, const T& sigma)
		{
			return _gearing * x 	* cva::exp(mu	* maturity) - _strike;
		}
	private:
		double _mu;
		double _gearing;
		double _strike;
		double _maturity;
	};

	template <typename T>
	T forwardFunction(const T& x,
		const T& mu, const T& sigma,
		const double a,
		const double b,
		const double maturity) 
	{
		return a * x 	* cva::exp(mu	* maturity ) - b;
	}

	//E[max(gearing * X(t)  - strike, 0) + shift] 
	class EuropeanCalculator {
	public:
		EuropeanCalculator(const double mu,
			const double gearing,
			const double strike,
			const double shift,
			const double maturity) : _mu(mu),
			_gearing(gearing), _strike(strike),
			_maturity(maturity), _shift(shift) {}
		template<typename T>
		T operator()(const T& x, const T& sigma) const
		{
			if (_maturity == 0.0) {
				return cva::zeroFloor(_gearing * x - _strike) + _shift;
			}
			T dplus = (cva::log(x / (_strike / _gearing))
				+ (_mu + sigma * sigma / 2.0) * _maturity)
				/ (sigma * std::sqrt(_maturity));
			T dminus = dplus
				- sigma * std::sqrt(_maturity);
			return (x * cva::exp(_mu * _maturity)
				* cva::normalCdf(dplus)
				+ cva::normalCdf(dminus) * (-_strike / _gearing)) 
				* _gearing + _shift;
		}
	private:
		double _mu;
		double _gearing;
		double _strike;
		double _maturity;
		double _shift;
	};

	template <typename T>
	T europeanFunction(const T& x,
		const T& mu, const T& sigma,
		const double a, const double b,
		const double c,
		const double maturity) 
	{
		if (maturity == 0.0) {
			return cva::zeroFloor(a * x - b) + c;
		}
		T dplus = (cva::log(x / (b / a))
			+ (mu + sigma * sigma / 2.0) * maturity) 
			/ (sigma * std::sqrt(maturity));
		T dminus = dplus 
			- sigma * std::sqrt(maturity);
		return (x * cva::exp(mu * maturity) 
			* cva::normalCdf(dplus)
			+ cva::normalCdf(dminus) * (-b / a)) * a + c;
	}

	//d/dx E[max(aX(t) + b, 0)]
	template <typename T>
	T europeanDelta(const T& x,
		const T& mu, const T& sigma,
		const double a, const double b,
		const double maturity)
	{
		T dplus = (cva::log(x / (b / a))
			+ (mu + sigma * sigma / 2.0) * maturity)
			/ (sigma * std::sqrt(maturity));
		T dminus = dplus - sigma 
			* std::sqrt(maturity);
		
		return cva::exp(mu * maturity)
			* cva::normalCdf(dplus) * a;
	}

	//E[max(aX(t) + b1, 0)] - E[max(-aX(t) - b2, 0)]
	template <typename T>
	T riskReversal(const T& x,
		const T& mu, const T& sigma,
		const double a, const double b1,
		const double b2,
		const double maturity)
	{
		//b1 > b2 >0
		return europeanFunction(x, mu, sigma, a, b1, 0.0,maturity)
			- europeanFunction(x, mu, sigma, -a, b2, 0.0, maturity);
	}

	//   E[max(aX(t) + b1, 0)] - E[max(aX(t) + b2, 0)]
	//+E[max(aX(t) + b1, 0)] - E[max(aX(t) + b2, 0)] + c
	template <typename T>
	T mountain(const T& x,
		const T& mu, const T& sigma,
		const double a, 
		const ublas::vector<double>& strikes,
		const double c,
		const double maturity)
	{
		//b1 < b2 < b3 < b4
		return europeanFunction(x, mu, sigma, a, strikes(0), 0.0, maturity)
			- europeanFunction(x, mu, sigma, a, strikes(1), 0.0, maturity)
			- europeanFunction(x, mu, sigma, a, strikes(2), 0.0, maturity)
			+ europeanFunction(x, mu, sigma, a, strikes(3), 0.0, maturity) + c;
	}
} // namespace cva