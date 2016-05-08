#include "stdafx.h"
#include "CppUnitTest.h"
#include <CvaCalculatorByAnalytic.h>
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	TEST_CLASS(Analytic)
	{
	public:
		TEST_METHOD(CvaAnalyticForward)
		{
			cva::Forward payoff(1.0, 100.0);
			double x0{ 100.0 };
			double sigma{ 0.3 };
			double mu{ 0.0 };

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const double dt = maturity / gridNum;

			cva::Dual<double> cva1 = cva::calcCvaFwdByAnalytic(x0, mu, sigma, payoff, maturity,
				gridNum, cva::undEnum);
			cva::Dual<double> cva2 = cva::calcCvaFwdByAnalytic(x0, mu, sigma, payoff, maturity,
				gridNum, cva::volEnum);

			Assert::AreEqual(
				// Expected value:
				226.18707825905656,
				// Actual value:
				cva1.value(),
				// Tolerance:
				1e-2,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				226.18707825905656,
				// Actual value:
				cva2.value(),
				// Tolerance:
				1e-2,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				6.6309353912952833,
				// Actual value:
				cva1.deriv(),
				// Tolerance:
				1e-2,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				722.205,
				// Actual value:
				cva2.deriv(),
				// Tolerance:
				1e-2,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}

		TEST_METHOD(CvaAnalyticEuropean)
		{
			cva::European payoff(1.0, 100.0);
			double x0{ 100.0 };
			double sigma{ 0.3 };
			double mu{ 0.0 };

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const double dt = maturity / gridNum;

			cva::Dual<double> cva1 = cva::calcCvaEurByAnalytic(x0, mu, sigma, payoff, maturity,
				gridNum, cva::undEnum);

			Assert::AreEqual(364.74370400275177, cva1.value(), 1e-2);
			Assert::AreEqual(6.8237185200137587, cva1.deriv(), 1e-2);

			cva::Dual<double> cva2 = cva::calcCvaEurByAnalytic(x0, mu, sigma, payoff, maturity,
				gridNum, cva::volEnum);

			Assert::AreEqual(364.74370400275177, cva2.value(), 1e-2);
			Assert::AreEqual(1127.3322640402168, cva2.deriv(), 1e-2);
		}

		TEST_METHOD(CvaAnalyticMountain)
		{
			cva::Mountain payoff(1.0, 90, 100, 110, 120);
			double x0{ 100.0 };
			double sigma{ 0.3 };
			double mu{ 0.0 };

			const double maturity = 10.0;
			const std::size_t gridNum = 10;
			const double dt = maturity / gridNum;

			cva::Dual<double> cva1 = cva::calcCvaMountainByAnalytic(x0, mu, sigma, payoff, maturity,
				gridNum, cva::undEnum);

			Assert::AreEqual(	7.0136937001802124, cva1.value(), 1e-2);
			Assert::AreEqual(0.038263062761274558, cva1.deriv(), 1e-2);

			cva::Dual<double> cva2 = cva::calcCvaMountainByAnalytic(x0, mu, sigma, payoff, maturity,
				gridNum, cva::volEnum);

			Assert::AreEqual(7.0136937001802124, cva2.value(), 1e-2);
			Assert::AreEqual(-28.496346740644185, cva2.deriv(), 1e-2);
		}
	};
}