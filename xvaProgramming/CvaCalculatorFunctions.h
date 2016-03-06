#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include "Regression.h"
#include "LsmFunctions.h"
#include <boost/bind.hpp>
#include "PathMaker.h"
#include "analytic_exposure_traits.h"
#include "CvaCalculator.h"
#include "RegressionExposure.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;

	//Calculate Cva By Analytic Exposure
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

	//Calculate Cva By Regression Exposure	
	template <typename T, typename U, typename P, typename C>
	T calcCvaByRegressionExposure(
		const Path<U>& pathForRegression,
		const Path<T>& pathForMonte,
		const PayOff<P>& payoff,
		const std::size_t numOfBasis,
		const CvaCalculator<C>& calculator)
	{
		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > lsmFunctions =
			makeLsmFunctions(numOfBasis, pathForRegression, payoff);
		Dual<double> cvaValue = calculator()(lsmFunctions, pathForMonte, payoff);
		return cvaValue;
	}


}//namespace cva
