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

	//Calculate Cva By Analytic Exposure
	//calclate df/dx, input x = Dual(x0, 1.0), mu = Dual(mu0, 0), sigma = Dual(sigma0, 0)
	//calclate df/dsigma, input x = Dual(x0, 0), mu = Dual(mu0, 0), sigma = Dual(sigma0, 1.0)
	template <typename T, typename P, typename C>
	T calcCvaByAnalyticExposure(
		const Path<T>& path,
		const PayOff<P>& payoff,
		const CvaCalculator<C>& calculator)
	{
		auto exposureFunctions
			= makeAnalyticExposureFunctions(path, payoff());
		T cvaValue = calculator()(
			exposureFunctions, path, payoff());
		return cvaValue;
	}

}//namespace cva
