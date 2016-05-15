#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include "Regression.h"
#include "lsm_exposure_traits.h"
#include <boost/bind.hpp>
#include "analytic_exposure_traits.h"
#include "CvaCalculator.h"
#include "Exposure.h"
#include "Regressor.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;

	//Calculate Cva By Analytic Exposure
	template <typename T, typename P, typename C>
	T calcCvaByAnalyticExposure(
		const Path<T>& path,
		const PayOff<P>& payoff,
		const CvaCalculator<C>& calculator)
	{
		AnalyticExposure<P> exposure(payoff());
		T cvaValue = calculator()(exposure, path, payoff());
		return cvaValue;
	}

	//Calculate Cva By Regression Exposure	
	template <typename T, typename R, typename S, typename B, typename P, typename C>
	T calcCvaByRegressionExposure(
		const Path<T>& path,
		const Regressor<R, P, S>& regressor,
		const ublas::vector<BasisFunctions<B>>& basisSeries,
		const PayOff<P>& payoff,
		const CvaCalculator<C>& calculator)
	{
		RegressionExposure<T, R, S, B, P> exposure(path, basisSeries, regressor);
		T cvaValue = calculator()(exposure, path, payoff());
		return cvaValue;
	}
}//namespace cva
