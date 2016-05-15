#include "stdafx.h"
#include "CppUnitTest.h"
#include <CvaCalculatorFunctions.h>
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	TEST_CLASS(AnalyticExposure)
	{
	public:
		TEST_METHOD(CvaAnalyticExposureForward)
		{
			cva::Forward payoff(1.0, 100.0);
			cva::Dual<double> x(100, 1.0);
			cva::Dual<double> sigma(0.3);
			cva::Dual<double> mu(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNum = 1000;
			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const cva::Path<cva::Dual<double>> path(
				x, mu, sigma, pathNum, gridNum, dt, seed);
		 
			cva::Dual<double> cva1 = cva::calcCvaByAnalyticExposure(
				path, payoff, cva::ExplicitCalculator());
			
			Assert::AreEqual(	236.96797947860327, cva1.value(), 1e-10);
			Assert::AreEqual(6.8166797947860260, cva1.deriv(), 1e-10);

			cva::Dual<double> cva2 = cva::calcCvaByAnalyticExposure(
				path, payoff, cva::ImplicitCalculator());

			Assert::AreEqual(212.437, cva2.value(), 1e-2);
			Assert::AreEqual(6.57137, cva2.deriv(), 1e-2);
		}

		TEST_METHOD(CvaAnalyticExposureEuropean)
		{
			cva::European payoff(1.0, 100.0);
			cva::Dual<double> x(100, 1.0);
			cva::Dual<double> sigma(0.3);
			cva::Dual<double> mu(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNum = 1000;
			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const cva::Path<cva::Dual<double>> path(
				x, mu, sigma, pathNum, gridNum, dt, seed);

			cva::Dual<double> cva1 = cva::calcCvaByAnalyticExposure(
				path, payoff, cva::ExplicitCalculator{});

			Assert::AreEqual(374.765, cva1.value(), 1e-2);
			Assert::AreEqual(6.98343, cva1.deriv(), 1e-2);

			cva::Dual<double> cva2 = cva::calcCvaByAnalyticExposure(
				path, payoff, cva::ImplicitCalculator{});

			Assert::AreEqual(350.204,	cva2.value(), 1e-2);
			Assert::AreEqual(6.64204, cva2.deriv(), 1e-2);
		}

		TEST_METHOD(CvaAnalyticExposureAsian)
		{
			cva::Asian payoff(1.0, 100.0);
			cva::Dual<double> x(100, 1.0);
			cva::Dual<double> sigma(0.3);
			cva::Dual<double> mu(0.0);

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const std::size_t pathNum = 1000;
			const double dt = maturity / gridNum;
			std::size_t seed = 1;
			const cva::Path<cva::Dual<double>> path(
				x, mu, sigma, pathNum, gridNum, dt, seed);

			cva::Dual<double> cva1 = cva::calcCvaByAnalyticExposure(
				path, payoff, cva::ExplicitCalculator{});

			Assert::AreEqual(172.15961465238428, cva1.value(), 1e-2);
			Assert::AreEqual(6.3285961465238341, cva1.deriv(), 1e-2);

			cva::Dual<double> cva2 = cva::calcCvaByAnalyticExposure(
				path, payoff, cva::ImplicitCalculator{});

			Assert::AreEqual(169.74907011287286, cva2.value(), 1e-2);
			Assert::AreEqual(6.3044907011287279, cva2.deriv(), 1e-2);
		}
	};
}