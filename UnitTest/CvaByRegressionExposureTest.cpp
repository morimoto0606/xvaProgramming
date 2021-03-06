#include "stdafx.h"
#include "CppUnitTest.h"
#include <Path.h>
#include <Payoff.h>
#include <CvaCalculator.h>
#include <CvaCalculatorFunctions.h>
#include <Regressor.h>
#include <BasisFunctions.h>
#include <boost/assign.hpp>
#include <fstream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace cva;
namespace UnitTest
{
	namespace ublas = boost::numeric::ublas;
	TEST_CLASS(CvaByLsm)
	{
	public:
		TEST_METHOD(CvaLsmExposureForwardCoeffFix)
		{
			Forward payoff(1.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 1000;
			const std::size_t pathNumForMonte = 1000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<double> pathForRegression(x0, mu0, sigma0, pathNumForRegression, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
						
			ublas::vector<BasisFunctions<double>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<double (double)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeriesForRegression(i) = cva::BasisFunctions<double>(basis);
			}
			Regressor<double, Forward> regressor(pathForRegression, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(214.46804145603323, cva1.value(), 1e-2); //value
				Assert::AreEqual(5.2087037857706706, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(214.46804145603323, cva2.value(), 1e-2); //value
				Assert::AreEqual(694.194, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(214.59570569475875, cva1.value(), 1e-2);
				Assert::AreEqual(5.4819570569475857, cva1.deriv(), 1e-2);
				Assert::AreEqual(214.59570569475875, cva2.value(), 1e-2); //value
				Assert::AreEqual(607.179, cva2.deriv(), 1e-2); //vega
			}
		}

		TEST_METHOD(CvaLsmExposureForwardCoeffShock)
		{
			Forward payoff(1.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 1000;
			const std::size_t pathNumForMonte = 1000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<Dual<double> > pathForRegression1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForRegression2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeriesForRegression(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}
			Regressor<Dual<double>, Forward> regressor1(pathForRegression1, basisSeriesForRegression, payoff);
			Regressor<Dual<double>, Forward> regressor2(pathForRegression2, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(214.46804145603323, cva1.value(), 1e-2); //value
				Assert::AreEqual(5.4806804145603829, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(214.46804145603323, cva2.value(), 1e-2); //value
				Assert::AreEqual(611.39430680391160, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(214.59570569475875, cva1.value(), 1e-2);
				Assert::AreEqual(5.4819570569475857, cva1.deriv(), 1e-2);
				Assert::AreEqual(214.59570569475875, cva2.value(), 1e-2); //value
				Assert::AreEqual(607.17858366809560, cva2.deriv(), 1e-2); //vega
			}
		}

		TEST_METHOD(CvaLsmExposureEuropeanCoeffFix)
		{
			European payoff(1.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 1000;
			const std::size_t pathNumForMonte = 1000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<double> pathForRegression(x0, mu0, sigma0, pathNumForRegression, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<BasisFunctions<double>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<double(double)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeriesForRegression(i) = cva::BasisFunctions<double>(basis);
			}
			Regressor<double, European> regressor(pathForRegression, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(362.59590176170042, cva1.value(), 1e-2); //value
				Assert::AreEqual(5.6457124778154464, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(362.59590176170042, cva2.value(), 1e-2); //value
				Assert::AreEqual(434.77677167369325, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(349.60846144728697, cva1.value(), 1e-2);
				Assert::AreEqual(6.6160846144728653, cva1.deriv(), 1e-2);
				Assert::AreEqual(349.60846144728697, cva2.value(), 1e-2); //value
				Assert::AreEqual(1017.1209824508499, cva2.deriv(), 1e-2); //vega
			}
		}

		TEST_METHOD(CvaLsmExposureEuropeanCoeffShock)
		{
			European payoff(1.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 1000;
			const std::size_t pathNumForMonte = 1000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<Dual<double> > pathForRegression1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForRegression2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeriesForRegression(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}
			Regressor<Dual<double>, European> regressor1(pathForRegression1, basisSeriesForRegression, payoff);
			Regressor<Dual<double>, European> regressor2(pathForRegression2, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(362.59590176170042, cva1.value(), 1e-2); //value
				Assert::AreEqual(6.7807220618627619, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(362.59590176170042, cva2.value(), 1e-2); //value
				Assert::AreEqual(1094.2243618553016, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(349.60846144728697, cva1.value(), 1e-2);
				Assert::AreEqual(6.6160846144728653, cva1.deriv(), 1e-2);
				Assert::AreEqual(349.60846144728697, cva2.value(), 1e-2); //value
				Assert::AreEqual(1017.1209824508499, cva2.deriv(), 1e-2); //vega
			}
		}

		TEST_METHOD(CvaLsmExposureAsianCoeffShockOneTimeBasis)
		{
			Asian payoff(1.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 10000;
			const std::size_t pathNumForMonte = 10000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<Dual<double> > pathForRegression1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForRegression2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeriesForRegression(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}
			Regressor<Dual<double>, Asian> regressor1(pathForRegression1, basisSeriesForRegression, payoff);
			Regressor<Dual<double>, Asian> regressor2(pathForRegression2, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(157.06369003327276, cva1.value(), 1e-2); //value
				Assert::AreEqual(6.0445369003326146, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(157.06369003327276, cva2.value(), 1e-2); //value
				Assert::AreEqual(839.87318628451112, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(158.16171646831756, cva1.value(), 1e-2); //value
				Assert::AreEqual(6.0555171646831756, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(158.16171646831756, cva2.value(), 1e-2); //value
				Assert::AreEqual(820.81286309983875, cva2.deriv(), 1e-2); //vega
			}
		}

		TEST_METHOD(CvaLsmExposureAsianCoeffShockAverageBasis)
		{
			typedef cva::Dual<double> state_type;
			typedef ublas::vector<state_type> basis_state_type;
			Asian payoff(1.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 10000;
			const std::size_t pathNumForMonte = 10000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<Dual<double> > pathForRegression1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForRegression2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<cva::BasisFunctions<basis_state_type>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<state_type(basis_state_type)>> basis(4);
				basis(0) = cva::PathwiseMonomial(i, 0);
				basis(1) = cva::PathwiseMonomial(i, 1);
				basis(2) = cva::PathwiseMonomial(i, 2);
				basis(3) = cva::TimewiseAverage(i, 1);
				basisSeriesForRegression(i) = cva::BasisFunctions<basis_state_type>(basis);
			}
			Regressor<state_type, Asian, basis_state_type> regressor1(pathForRegression1, basisSeriesForRegression, payoff);
			Regressor<state_type, Asian, basis_state_type> regressor2(pathForRegression2, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<basis_state_type>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<state_type(basis_state_type)>> basis(4);
				basis(0) = cva::PathwiseMonomial(i, 0);
				basis(1) = cva::PathwiseMonomial(i, 1);
				basis(2) = cva::PathwiseMonomial(i, 2);
				basis(3) = cva::TimewiseAverage(i, 1);
				basisSeries(i) = cva::BasisFunctions<basis_state_type>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(163.30274506468109, cva1.value(), 1e-2); //value
				Assert::AreEqual(6.2127621326887192, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(163.30274506468109, cva2.value(), 1e-2); //value
				Assert::AreEqual(809.61260277154497, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(166.78122154356896, cva1.value(), 1e-2); //value
				Assert::AreEqual(6.4038122154356953, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(166.78122154356896, cva2.value(), 1e-2); //value
				Assert::AreEqual(830.62042172184510, cva2.deriv(), 1e-2); //vega
			}
		}
		TEST_METHOD(CvaLsmExposureTrfCoeffShockOneTimeBasis)
		{
			Trf payoff(1.0, 0.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 10000;
			const std::size_t pathNumForMonte = 10000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<Dual<double> > pathForRegression1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForRegression2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeriesForRegression(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}
			Regressor<Dual<double>, Trf> regressor1(pathForRegression1, basisSeriesForRegression, payoff);
			Regressor<Dual<double>, Trf> regressor2(pathForRegression2, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<cva::Dual<double>>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<cva::Dual<double>(cva::Dual<double>)>> basis(3);
				basis(0) = cva::Monomial(0);
				basis(1) = cva::Monomial(1);
				basis(2) = cva::Monomial(2);
				basisSeries(i) = cva::BasisFunctions<cva::Dual<double>>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(702.37882842427825, cva1.value(), 1e-2); //value
				//Assert::AreEqual(6.2346689222621459, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(702.37882842427825, cva2.value(), 1e-2); //value
				//Assert::AreEqual(861.79994612248208, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(680.19411364190103, cva1.value(), 1e-2); //value
				//Assert::AreEqual(6.2407073612623662, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(680.19411364190103, cva2.value(), 1e-2); //value
				//Assert::AreEqual(838.63513538016389, cva2.deriv(), 1e-2); //vega
			}
		}

		TEST_METHOD(CvaLsmExposureTrfCoeffShockOneTimePathwiseBasis)
		{
			typedef cva::Dual<double> state_type;
			typedef ublas::vector<state_type> basis_state_type;
			Trf payoff(1.0, 0.0, 100.0);
			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 1000;
			const std::size_t pathNumForMonte = 1000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<Dual<double> > pathForRegression1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForRegression2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<cva::BasisFunctions<basis_state_type>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<state_type(basis_state_type)>> basis(5);
				basis(0) = cva::PathwiseMonomial(i, 0);
				basis(1) = cva::PathwiseMonomial(i, 1);
				basis(2) = cva::PathwiseMonomial(i, 2);
				basis(3) = cva::TimewiseAverage(i, 1);
				basis(4) = cva::TimewiseAverage(i, 2);
				basisSeriesForRegression(i) = cva::BasisFunctions<basis_state_type>(basis);
			}
			Regressor<state_type, Trf, basis_state_type> regressor1(pathForRegression1, basisSeriesForRegression, payoff);
			Regressor<state_type, Trf, basis_state_type> regressor2(pathForRegression2, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<basis_state_type>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<state_type(basis_state_type)>> basis(5);
				basis(0) = cva::PathwiseMonomial(i, 0);
				basis(1) = cva::PathwiseMonomial(i, 1);
				basis(2) = cva::PathwiseMonomial(i, 2);
				basis(3) = cva::TimewiseAverage(i, 1);
				basis(4) = cva::TimewiseAverage(i, 2);
				basisSeries(i) = cva::BasisFunctions<basis_state_type>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ExplicitCalculator{});

				Assert::AreEqual(702.37882842427825, cva1.value(), 1e-2); //value
																		  //Assert::AreEqual(6.2346689222621459, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(702.37882842427825, cva2.value(), 1e-2); //value
																		  //Assert::AreEqual(861.79994612248208, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(680.19411364190103, cva1.value(), 1e-2); //value
																		  //Assert::AreEqual(6.2407073612623662, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(680.19411364190103, cva2.value(), 1e-2); //value
																		  //Assert::AreEqual(838.63513538016389, cva2.deriv(), 1e-2); //vega
			}
		}
		TEST_METHOD(CvaLsmExposureTrfCoeffShockAverageBasis)
		{
			typedef cva::Dual<double> state_type;
			typedef ublas::vector<state_type> basis_state_type;
			Trf payoff(1.0, 0.0, 100.0);

			double x0 = 100.0;
			cva::Dual<double> x1(100.0);
			cva::Dual<double> x2(100.0, 1.0);
			double sigma0 = 0.3;
			cva::Dual<double> sigma1(0.3);
			cva::Dual<double> sigma2(0.3, 1.0);
			double mu0 = 0.0;
			cva::Dual<double> mu1(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNumForRegression = 1000;
			const std::size_t pathNumForMonte = 1000;

			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const Path<Dual<double> > pathForRegression1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForRegression2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte1(x2, mu1, sigma1, pathNumForMonte, gridNum, dt, seed);
			const Path<Dual<double> > pathForMonte2(x1, mu1, sigma2, pathNumForMonte, gridNum, dt, seed);

			ublas::vector<cva::BasisFunctions<basis_state_type>> basisSeriesForRegression(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<state_type(basis_state_type)>> basis(3);
				basis(0) = cva::TimewiseAverage(i, 0);
				basis(1) = cva::TimewiseAverage(i, 1);
				basis(2) = cva::TimewiseAverage(i, 2);
				basisSeriesForRegression(i) = cva::BasisFunctions<basis_state_type>(basis);
			}
			Regressor<state_type, Trf, basis_state_type> regressor1(pathForRegression1, basisSeriesForRegression, payoff);
			Regressor<state_type, Trf, basis_state_type> regressor2(pathForRegression2, basisSeriesForRegression, payoff);

			ublas::vector<cva::BasisFunctions<basis_state_type>> basisSeries(gridNum);
			for (std::size_t i = 0; i < gridNum; ++i) {
				ublas::vector<boost::function<state_type(basis_state_type)>> basis(3);
				basis(0) = cva::TimewiseAverage(i, 0);
				basis(1) = cva::TimewiseAverage(i, 1);
				basis(2) = cva::TimewiseAverage(i, 2);
				basisSeries(i) = cva::BasisFunctions<basis_state_type>(basis);
			}

			//Explicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ExplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ExplicitCalculator{});

				//Assert::AreEqual(705.52326629895344, cva1.value(), 1e-2); //value
				//Assert::AreEqual(6.3997144927420448, cva1.deriv(), 1e-2); //delta 
				//Assert::AreEqual(705.52326629895344, cva2.value(), 1e-2); //value
				//Assert::AreEqual(840.37419651279436, cva2.deriv(), 1e-2); //vega
			}
			//Implicit
			{
				Dual<double> cva1 = calcCvaByRegressionExposure(
					pathForMonte1, regressor1, basisSeries, payoff, cva::ImplicitCalculator{});
				Dual<double> cva2 = calcCvaByRegressionExposure(
					pathForMonte2, regressor2, basisSeries, payoff, cva::ImplicitCalculator{});

				Assert::AreEqual(680.30549423151172, cva1.value(), 1e-2); //value
				//::AreEqual(680.30549423151172, cva1.deriv(), 1e-2); //delta 
				Assert::AreEqual(680.30549423151172, cva2.value(), 1e-2); //value
				//Assert::AreEqual(830.62042172184510, cva2.deriv(), 1e-2); //vega
			}
		}
	};
}