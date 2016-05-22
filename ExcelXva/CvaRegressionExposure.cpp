#include "stdafx.h"
#include <windows.h>
#undef min
#undef max
#include "CvaRegressionExposure.h"

#include <Path.h>
#include <Payoff.h>
#include <CvaCalculator.h>
#include <CvaCalculatorFunctions.h>
#include <Regressor.h>
#include <BasisFunctions.h>
#include <boost/assign.hpp>

using namespace cva;
double __stdcall  calcCvaLsmExposureEuropean(
	const double spot,
	const double strike,
	const double vol,
	const double maturity,
	const int gridNum,
	const long seed,
	const long pathNumForRegression,
	const long pathNumForMonte,
	const int orderOfBasisFunction,
	const bool isExplicit,
	const bool isCoeffShock)
{
	//,
		//const long seed,
		//const long pathNumForRegression,
		//const long pathNumForMonte,
		//const int orderOfBasisFunction,
		//const bool isExplicit,
		//const bool isCoeffShock,
		//const int calcType


		const int calcType = 0;

	namespace ublas = boost::numeric::ublas;
	European payoff{ 1.0, strike };

	cva::Dual<double> x = calcType != 2
		? cva::Dual<double>(spot, 1.0) 
		: cva::Dual<double>(spot, 0.0);
	cva::Dual<double> sigma = calcType != 2
		? cva::Dual<double>(vol, 0.0)
		: cva::Dual<double>(vol, 1.0);

	cva::Dual<double> mu(0.0);
	
	const double dt = maturity / gridNum;
	double result = 0.0;
	
	if (!isCoeffShock) {
		const Path<double> pathForRegression(
			x.value(), mu.value(), sigma.value(), pathNumForRegression, gridNum, dt, seed);

		ublas::vector<BasisFunctions<double>> basisSeriesForRegression(gridNum);
		for (std::size_t i = 0; i < gridNum; ++i) {
			ublas::vector<boost::function<double(double)>> basis(orderOfBasisFunction + 1);
			for (int j = 0; j <= orderOfBasisFunction; ++j) {
				basis(j) = cva::Monomial(j);
			}
			basisSeriesForRegression(i) = cva::BasisFunctions<double>(basis);
		}
		Regressor<double, European> regressor(pathForRegression, basisSeriesForRegression, payoff);

		const Path<Dual<double> > pathForMonte(
			x, mu, sigma, pathNumForMonte, gridNum, dt, seed);

		ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
		for (std::size_t i = 0; i < gridNum; ++i) {
			ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>>
				basis(orderOfBasisFunction + 1);
			for (int j = 0; j <= orderOfBasisFunction; ++j) {
				basis(j) = cva::Monomial(j);
			}
			basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
		}

		Dual<double> cva = isExplicit
			? calcCvaByRegressionExposure(
				pathForMonte, regressor, basisSeries, payoff, cva::ExplicitCalculator{})
			: calcCvaByRegressionExposure(
				pathForMonte, regressor, basisSeries, payoff, cva::ImplicitCalculator{});

		result = (calcType == 0) ? cva.value() : cva.deriv();
	}
	else {
		const Path<Dual<double> > pathForRegression(
			x, mu, sigma, pathNumForMonte, gridNum, dt, seed);

		ublas::vector<cva::BasisFunctions<cva::Dual<double>>>
			basisSeriesForRegression(gridNum);
		for (std::size_t i = 0; i < gridNum; ++i) {
			ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> 	
				basis(orderOfBasisFunction + 1);
			for (int j = 0; j <= orderOfBasisFunction; ++j) {
				basis(j) = cva::Monomial(j);
			}
			basisSeriesForRegression(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
		}
		Regressor<Dual<double>, European> regressor(
			pathForRegression, basisSeriesForRegression, payoff);

		const Path<Dual<double> > pathForMonte(
			x, mu, sigma, pathNumForMonte, gridNum, dt, seed);

		ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
		for (std::size_t i = 0; i < gridNum; ++i) {
			ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> 			
				basis(orderOfBasisFunction + 1);
			for (int j = 0; j <= orderOfBasisFunction; ++j) {
				basis(j) = cva::Monomial(j);
			}
			basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
		}

		Dual<double> cva = isExplicit
			? calcCvaByRegressionExposure(
				pathForMonte, regressor, basisSeries, payoff, cva::ExplicitCalculator{})
			: calcCvaByRegressionExposure(
				pathForMonte, regressor, basisSeries, payoff, cva::ImplicitCalculator{});

		result = (calcType == 0) ? cva.value() : cva.deriv();
	}
	return result;
}