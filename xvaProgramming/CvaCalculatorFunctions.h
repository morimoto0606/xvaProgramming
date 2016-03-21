#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include "Regression.h"
#include "lsm_exposure_traits.h"
#include <boost/bind.hpp>
#include "PathMaker.h"
#include "analytic_exposure_traits.h"
#include "CvaCalculator.h"
#include "Exposure.h"

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
	template <typename T, typename U, typename D, 
		typename R, typename P, typename C>
	T calcCvaByRegressionExposure(
		const Path<U>& pathForRegression,
		const Path<T>& pathForMonte,
		const ublas::vector<BasisFunctions<D, R>>& basisSeries,
		const PayOff<P>& payoff,
		const CvaCalculator<C>& calculator)
	{
		RegressionExposure<T, D, R, P> exposure(
			pathForRegression, basisSeries, payoff());
		T cvaValue = calculator()(exposure, pathForMonte, payoff());
		return cvaValue;
	}


}//namespace cva
