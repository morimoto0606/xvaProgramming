#include "stdafx.h"
#include "CppUnitTest.h"
#include <AnalyticExposure.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	TEST_CLASS(AnalyticExposure)
	{
	public:
		TEST_METHOD(AnalyticExposureForward)
		{
			//cva::Forward payoff(1.0, 100.0);
			//cva::Dual<double> cvaForward = cva::calcCvaByAnalyticExposure(
			//	100, 0.0, 0.3, payoff, 10.0, 10, 1000, 1,
			//	productTypeEnum::eurEnum, false);
			//Assert::AreEqual(
			//	// Expected value:
			//	x1 + x2,
			//	// Actual value:
			//	(dual1 + dual2).value(),
			//	// Tolerance:
			//	1e-10,
			//	// Message:
			//	L"Basic test failed",
			//	// Line number - used if there is no PDB file:
			//	LINE_INFO());

			//Assert::AreEqual(
			//	// Expected value:
			//	y1 + y2,
			//	// Actual value:
			//	(dual1 + dual2).deriv(),
			//	// Tolerance:
			//	1e-10,
			//	// Message:
			//	L"Basic test failed",
			//	// Line number - used if there is no PDB file:
			//	LINE_INFO());
		}

		TEST_METHOD(DualMinus)
		{
			double x1 = 1.1543543543;
			double x2 = 2.14535435;
			double y1 = -2.456;
			double y2 = 5.03343;
			cva::Dual<double> dual1(x1, y1);
			cva::Dual<double> dual2(x2, y2);

			cva::Dual<double> minus = dual1 - dual2;

			Assert::AreEqual(
				// Expected value:
				x1 - x2,
				// Actual value:
				(dual1 - dual2).value(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				y1 - y2,
				// Actual value:
				(dual1 - dual2).deriv(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}

		TEST_METHOD(DualProd)
		{
			double x1 = 1.1543543543;
			double x2 = 2.14535435;
			double y1 = -2.456;
			double y2 = 5.03343;
			cva::Dual<double> dual1(x1, y1);
			cva::Dual<double> dual2(x2, y2);

			cva::Dual<double> prod = dual1 * dual2;

			Assert::AreEqual(
				// Expected value:
				x1 * x2,
				// Actual value:
				prod.value(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				y1 * x2 + x1 * y2,
				// Actual value:
				prod.deriv(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}

		TEST_METHOD(DivProd)
		{
			double x1 = 1.1543543543;
			double x2 = 2.14535435;
			double y1 = -2.456;
			double y2 = 5.03343;
			cva::Dual<double> dual1(x1, y1);
			cva::Dual<double> dual2(x2, y2);

			cva::Dual<double> div = dual1 / dual2;

			Assert::AreEqual(
				// Expected value:
				x1 / x2,
				// Actual value:
				div.value(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				(y1 * x2 - x1 * y2) / (x2 * x2),
				// Actual value:
				div.deriv(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}
	};
}