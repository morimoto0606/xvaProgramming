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
		TEST_METHOD(AnalyticExposureForward)
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
			cva::ExplicitCalculator calculator;
			cva::Dual<double> cvaForward = cva::calcCvaByAnalyticExposure(
				path, payoff, calculator);
			std::cout << cvaForward.value() << ',' << cvaForward.deriv() << std::endl;
			
			Assert::AreEqual(
				// Expected value:
				236.96797947860327,
				// Actual value:
				cvaForward.value(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				6.8166797947860260,
				// Actual value:
				cvaForward.deriv(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}
	};
}