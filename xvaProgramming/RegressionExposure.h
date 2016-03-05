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

namespace cva {
	namespace ublas = boost::numeric::ublas;
	//Calculate Cva By Regression ExposureWithModification
	template <typename P, typename C, typename O>
	Dual<double> calcCvaByRegressionExposure(
		const double x0, const double mu, const double sigma,
		const PayOff<P>& payoff, const double maturity,
		const std::size_t gridNum, const std::size_t pathNumForRegression,
		const std::size_t pathNumForMonte,
		const std::size_t numOfBasis,
		const shockTypeEnum shockType, const std::size_t seed,
		const bool useImplicitMethod, const bool isCoeffShock,
		const O& analyticOptionFunction)
	{
		const double dt = maturity / gridNum;
		const Path<C> pathForRegression
			= makePath<C>(x0, mu, sigma, dt, gridNum,
				pathNumForRegression, shockType, seed);

		const Path<Dual<double> > pathForMonte
			= makePath<Dual<double> >(x0, mu, sigma, dt, gridNum,
				pathNumForMonte, shockType, seed);
		double trueValue = analyticOptionFunction(x0, sigma);
		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > lsmFunctions =
			makeLsmFunctions(numOfBasis, 
				pathForRegression, payoff, trueValue);

		Dual<double> cvaValue = useImplicitMethod
			? calcCvaUsingImplicitExposure(lsmFunctions, pathForMonte, payoff, dt)
			: calcCvaUsingExplicitExposure(lsmFunctions, pathForMonte, dt);
		return cvaValue;
	}

	//Calculate Cva By Regression Exposure	
	template <typename P, typename C>
	Dual<double> calcCvaByRegressionExposure(
		const double x0, const double mu, const double sigma,
		const PayOff<P>& payoff, const double maturity,
		const std::size_t gridNum, const std::size_t pathNumForRegression,
		const std::size_t pathNumForMonte,
		const std::size_t numOfBasis,
		const shockTypeEnum shockType, const std::size_t seed,
		const bool useImplicitMethod, const bool isCoeffShock)
	{
		const double dt = maturity / gridNum;
		const Path<C> pathForRegression
			= makePath<C>(x0, mu, sigma, dt, gridNum,
				pathNumForRegression, shockType, seed);

		const Path<Dual<double> > pathForMonte
			= makePath<Dual<double> >(x0, mu, sigma, dt, gridNum,
				pathNumForMonte, shockType, seed);

		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > lsmFunctions =
			makeLsmFunctions(numOfBasis, pathForRegression, payoff);

		Dual<double> cvaValue = useImplicitMethod
			? calcCvaUsingImplicitExposure(lsmFunctions, pathForMonte, payoff, dt)
			: calcCvaUsingExplicitExposure(lsmFunctions, pathForMonte, dt);
		return cvaValue;
	}

	template <typename P, typename C>
	ublas::vector<boost::function<Dual<double>(
		const Dual<double>&)> > makeLsmFunctions(
			std::size_t numOfBasis,
			const Path<C>& path, const PayOff<P>& payoff,
			const double trueValueForModification)
	{
		// lsmFunctions(i),  i = 0, 1, ...,  gridNum 
		ublas::vector<boost::function<C(const C&)> >
			functions(numOfBasis);
		for (std::size_t i = 0; i < numOfBasis; ++i) {
			functions(i) = boost::function<C(const C&)>
				(Monomial(static_cast<double>(i)));
		}
		ublas::vector<boost::function<
			Dual<double>(const Dual<double>&)> >  dualFunctions(numOfBasis);
		for (std::size_t i = 0; i < numOfBasis; ++i) {
			dualFunctions(i) = boost::function<
				Dual<double>(const Dual<double>&)>
				(Monomial(static_cast<double>(i)));
		}

		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> >
			lsmFunctions(path.gridNum() + 1);
		for (std::size_t gridIndex = 0; gridIndex <= path.gridNum(); ++gridIndex) {
			ublas::vector<C> coeff;
			coeff = regresssion(gridIndex, payoff, path, functions);
			Dual<double> v0 = 
				LsmFunction<Dual<double>, C>(
					coeff, dualFunctions)(Dual<double>(path.getPathValue(0, 0)));
			const double modificationValue 
				= trueValueForModification - v0.value();
			coeff(0) = coeff(0) + modificationValue;
			lsmFunctions(gridIndex)
				= LsmFunction<Dual<double>, C>(
					coeff, dualFunctions);
		}
		return lsmFunctions;
	}

	template <typename P, typename C>
	ublas::vector<boost::function<Dual<double>(
		const Dual<double>&)> > makeLsmFunctions(
			const std::size_t numOfBasis,
			const Path<C>& path, const PayOff<P>& payoff)
	{
		return makeLsmFunctions(numOfBasis, path, payoff, 0.0);
	}
} //namespace cva {
