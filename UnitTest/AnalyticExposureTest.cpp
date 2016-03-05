#include "stdafx.h"
#include "CppUnitTest.h"
#include <AnalyticExposure.h>
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
			cva::ExplicitCalculator calculator;
			cva::Dual<double> cvaForward = cva::calcCvaByAnalyticExposure(
				x, mu, sigma, payoff, 10.0, 10, 1000, 1, calculator);
			std::cout << cvaForward.value() << ',' << cvaForward.deriv() << std::endl;
			
			Assert::AreEqual(
				// Expected value:
				271.98840155308841,
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
				6.4808840155308722,
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