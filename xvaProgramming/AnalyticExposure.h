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

namespace cva {
	namespace ublas = boost::numeric::ublas;


	//Calculate Cva By Analytic Exposure
	template <typename T>
	Dual<double> calcCvaByAnalyticExposure(
		const double x0, const double mu, const double sigma,
		const PayOff<T>& payoff, const double maturity,
		const std::size_t gridNum, const std::size_t pathNum,
		const std::size_t seed,
		const bool useImplicitMethod)
	{
		const double dt = maturity / gridNum;

		const Path<Dual<double> > pathDual = makePath<Dual<double> >(
			x0, mu, sigma, dt, gridNum, pathNum, shockType, seed);

		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > exposureFunctions
			= analytic_exposure_traits<PayOff<T>>::apply(
			mu, sigma, pathDual.gridNum(), maturity, payoff);
		
		Dual<double> cvaValue = useImplicitMethod
			? calcCvaUsingImplicitExposure(
				exposureFunctions, pathDual, payoff, dt)
			: calcCvaUsingExplicitExposure(
				exposureFunctions, pathDual, dt);
		return cvaValue;
	}

	ublas::vector<boost::function<Dual<double>(
		const Dual<double>&)> >
		makeAnalyticForwardExposures(
			const double mu, const double sigma,
			const std::size_t gridNum, const double maturity,
			const double fwdCoeffA, const double fwdCoeffB)
	{
		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > fwdFunctions(gridNum + 1);
		for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
			const double tau = maturity
				* static_cast<double>(gridNum - gridIndex)
				/ static_cast<double>(gridNum);
			fwdFunctions(gridIndex)
				= boost::bind(forwardFunction<Dual<double> >, _1,
					Dual<double>(mu), Dual<double>(sigma),
					fwdCoeffA, fwdCoeffB, tau);
		}
		return fwdFunctions;
	}

	ublas::vector<boost::function<Dual<double>(
		const Dual<double>&)> >
		makeAnalyticEuropeanExposures(
			const double mu, const double sigma,
			const std::size_t gridNum, const double maturity,
			const double eurCoeffA, const double eurCoeffB,
			const double eurCoeffC)
	{
		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > eurFunctions(gridNum + 1);
		for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
			const double tau = maturity
				* static_cast<double>(gridNum - gridIndex)
				/ static_cast<double>(gridNum);
			eurFunctions(gridIndex)
				= boost::bind(europeanFunction<Dual<double> >, _1,
					Dual<double>(mu), Dual<double>(sigma),
					eurCoeffA, eurCoeffB, eurCoeffC, tau);
		}
		return eurFunctions;
	}

	ublas::vector<boost::function<Dual<double>(
		const Dual<double>&)> >
		makeAnalyticMountainExposures(
			const double mu, const double sigma,
			const std::size_t gridNum, const double maturity,
			const double gearing, 
			const ublas::vector<double>& strikes,
			const double payoffShift)
	{
		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > functions(gridNum + 1);
		for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
			const double tau = maturity
				* static_cast<double>(gridNum - gridIndex)
				/ static_cast<double>(gridNum);
			functions(gridIndex)
				= boost::bind(mountain<Dual<double> >, _1,
					Dual<double>(mu), Dual<double>(sigma),
					gearing, strikes, 
					payoffShift, tau);
		}
		return functions;
	}

	ublas::vector<boost::function<Dual<double>(
		const Dual<double>&)> >
		makeAnalyticRiskReversalExposures(
			const double mu, const double sigma,
			const std::size_t gridNum, const double maturity,
			const double gearing, const double strike1,
			const double strike2)
	{
		ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > rrFunctions(gridNum);
		for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
			const double tau = maturity
				* static_cast<double>(gridNum - gridIndex)
				/ static_cast<double>(gridNum);
			rrFunctions(gridIndex)
				= boost::bind(riskReversal<Dual<double> >, _1,
					Dual<double>(mu), Dual<double>(sigma),
					gearing, strike1, strike2, tau);
		}
		return rrFunctions;
	}
}//namespace cva
