#include "stdafx.h"
#include "CppUnitTest.h"
#include <Path.h>
#include <Payoff.h>
#include <CvaCalculator.h>
#include <CvaCalculatorFunctions.h>
#include <BasisFunctions.h>
#include <boost/assign.hpp>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace cva;
namespace UnitTest
{
	namespace ublas = boost::numeric::ublas;
	TEST_CLASS(CvaByLsm)
	{
	public:
		TEST_METHOD(CvaLsmForwardExplicit)
		{
			Forward payoff(1.0, 100.0);
			double x0 = 100;
			double sigma0 = 0.3;
			double mu0 = 0.0;
			cva::Dual<double> x(x0, 1.0);
			cva::Dual<double> sigma(sigma0);
			cva::Dual<double> mu(mu0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 1000;
			const std::size_t pathNumForMonte = 1000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<double> pathForRegression(x0, mu0, sigma0, pathNumForRegression, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte(x, mu, sigma, pathNumForMonte, gridNum, dt, seed);
			cva::ExplicitCalculator calculator;
			
			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}
			Dual<double> cva = calcCvaByRegressionExposure(
				pathForMonte, pathForMonte, basisSeries, payoff, calculator);

			Assert::AreEqual(214.46804145603323, cva.value(), 1e-2);
			Assert::AreEqual(5.4806804145603829, cva.deriv(), 1e-2);
		}

		//TEST_METHOD(CvaLsmForwardImplicit)
		//{
		//	Forward payoff(1.0, 100.0);
		//	double x0 = 100;
		//	double sigma0 = 0.3;
		//	double mu0 = 0.0;
		//	cva::Dual<double> x(x0, 1.0);
		//	cva::Dual<double> sigma(sigma0);
		//	cva::Dual<double> mu(mu0);

		//	const double maturity = 10.0;
		//	const std::size_t gridNum = 10;
		//	const std::size_t pathNumForRegression = 1000;
		//	const std::size_t pathNumForMonte = 1000;

		//	const double dt = maturity / gridNum;
		//	std::size_t seed = 1;
		//	const Path<double> pathForRegression(x0, mu0, sigma0, pathNumForRegression, gridNum, dt, seed);
		//	const Path<Dual<double> > pathForMonte(x, mu, sigma, pathNumForMonte, gridNum, dt, seed);
		//	cva::ImplicitCalculator calculator;
		//	ublas::vector<cva::BasisFunctions<cva::Dual<double>, cva::Dual<double>>> basisSeries(gridNum);
		//	for (std::size_t i = 0; i < gridNum; ++i) {
		//		ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
		//		basis(0) = cva::Monomial(0);
		//		basis(1) = cva::Monomial(1);
		//		basis(2) = cva::Monomial(2);
		//		basisSeries(i) = cva::BasisFunctions<cva::Dual<double>, cva::Dual<double>>(basis);
		//	}
		//	Dual<double> cva = calcCvaByRegressionExposure(
		//		pathForRegression, pathForMonte, basisSeries, payoff, calculator);

		//	Assert::AreEqual(248.904, cva.value(), 1e-2);
		//	Assert::AreEqual(6.25004, cva.deriv(), 1e-2);
		//}

		//TEST_METHOD(CvaLsmEuropeanExplicit)
		//{
		//	European payoff(1.0, 100.0);
		//	double x0 = 100;
		//	double sigma0 = 0.3;
		//	double mu0 = 0.0;
		//	cva::Dual<double> x(x0, 1.0);
		//	cva::Dual<double> sigma(sigma0);
		//	cva::Dual<double> mu(mu0);

		//	const double maturity = 10.0;
		//	const std::size_t gridNum = 10;
		//	const std::size_t pathNumForRegression = 1000;
		//	const std::size_t pathNumForMonte = 1000;

		//	const double dt = maturity / gridNum;
		//	std::size_t seed = 1;
		//	const Path<double> pathForRegression(x0, mu0, sigma0, pathNumForRegression, gridNum, dt, seed);
		//	const Path<Dual<double> > pathForMonte(x, mu, sigma, pathNumForMonte, gridNum, dt, seed);
		//	cva::ExplicitCalculator calculator;
		//	ublas::vector<cva::BasisFunctions<cva::Dual<double>, cva::Dual<double>>> basisSeries(gridNum);
		//	for (std::size_t i = 0; i < gridNum; ++i) {
		//		ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
		//		basis(0) = cva::Monomial(0);
		//		basis(1) = cva::Monomial(1);
		//		basis(2) = cva::Monomial(2);
		//		basisSeries(i) = cva::BasisFunctions<cva::Dual<double>, cva::Dual<double>>(basis);
		//	}
		//	Dual<double> cva = calcCvaByRegressionExposure(pathForRegression, pathForMonte, basisSeries, payoff, calculator);

		//	Assert::AreEqual(189.249, cva.value(), 1e-2);
		//	Assert::AreEqual(4.69081, cva.deriv(), 1e-2);
		//}

		//TEST_METHOD(CvaLsmEuropeanImplicit)
		//{
		//	European payoff(1.0, 100.0);
		//	double x0 = 100;
		//	double sigma0 = 0.3;
		//	double mu0 = 0.0;
		//	cva::Dual<double> x(x0, 1.0);
		//	cva::Dual<double> sigma(sigma0);
		//	cva::Dual<double> mu(mu0);

		//	const double maturity = 10.0;
		//	const std::size_t gridNum = 10;
		//	const std::size_t pathNumForRegression = 1000;
		//	const std::size_t pathNumForMonte = 1000;

		//	const double dt = maturity / gridNum;
		//	std::size_t seed = 1;
		//	const Path<double> pathForRegression(x0, mu0, sigma0, pathNumForRegression, gridNum, dt, seed);
		//	const Path<Dual<double> > pathForMonte(x, mu, sigma, pathNumForMonte, gridNum, dt, seed);
		//	cva::ImplicitCalculator calculator;
		//	ublas::vector<cva::BasisFunctions<cva::Dual<double>, cva::Dual<double>>> basisSeries(gridNum);
		//	for (std::size_t i = 0; i < gridNum; ++i) {
		//		ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
		//		basis(0) = cva::Monomial(0);
		//		basis(1) = cva::Monomial(1);
		//		basis(2) = cva::Monomial(2);
		//		basisSeries(i) = cva::BasisFunctions<cva::Dual<double>, cva::Dual<double>>(basis);
		//	}

		//	Dual<double> cva = calcCvaByRegressionExposure(
		//		pathForRegression, pathForMonte, basisSeries, payoff, calculator);

		//	Assert::AreEqual(300.122, cva.value(), 1e-2);
		//	Assert::AreEqual(5.42622, cva.deriv(), 1e-2);
		//}

	};
}