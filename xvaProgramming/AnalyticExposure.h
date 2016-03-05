#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include "Regression.h"
#include "LsmFunction.h"
#include <boost/bind.hpp>
#include "PathMaker.h"
#include "analytic_exposure_traits.h"
#include "CvaCalculator.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;

	template <typename T, typename P>
	ublas::vector<boost::function<T (const T&)>> makeAnalyticExposureFunctions
		(const P& payoff, const const T& mu, const T& sigma, 
			const std::size_t gridNum, const double maturity)
	{
		return analytic_exposure_traits<P>::apply(
			mu, sigma, gridNum, maturity, payoff);
	}

	//Calculate Cva By Analytic Exposure
	//calclate df/dx, input x = Dual(x0, 1.0), mu = Dual(mu0, 0), sigma = Dual(sigma0, 0)
	//calclate df/dsigma, input x = Dual(x0, 0), mu = Dual(mu0, 0), sigma = Dual(sigma0, 1.0)
	template <typename T, typename P, typename C>
	T calcCvaByAnalyticExposure(
		const T& x, const T& mu, const T& sigma,
		const P& payoff, const double maturity,
		const std::size_t gridNum, const std::size_t pathNum,
		const std::size_t seed,
		const CvaCalculator<C>& calculator)
	{
		const double dt = maturity / gridNum;
		const Path<T> path(x, mu, sigma, pathNum, gridNum, dt, seed);

		ublas::vector<boost::function<T(
			const T&)> > exposureFunctions
			= makeAnalyticExposureFunctions(
				payoff, mu, sigma, gridNum, maturity);
		
		T cvaValue = calculator()(
			exposureFunctions, path, payoff);
		return cvaValue;
	}

}//namespace cva
