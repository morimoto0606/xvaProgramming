#include "stdafx.h"
#include "CppUnitTest.h"
#include <function_traits.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	TEST_CLASS(function_traits)
	{
	public:
		TEST_METHOD(zero_floor_traits)
		{
			double x = 1.1543543543;
			double y = -2.456;
			cva::Dual<double> dual(x, y);

			double x2 = -1.1543543543;
			double y2 = -2.456;
			cva::Dual<double> dual2(x2, y2);

			Assert::AreEqual(
				// Expected value:
				std::max(x, 0.0),
				// Actual value:
				cva::zero_floor_traits<double>().apply(x),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				std::max(x, 0.0),
				// Actual value:
				cva::zero_floor_traits<cva::Dual<double>>().apply(dual).value(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		
			Assert::AreEqual(
				// Expected value:
				x > 0 ? y : 0.0,
				// Actual value:
				cva::zero_floor_traits<cva::Dual<double>>().apply(dual).deriv(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				std::max(x2, 0.0),
				// Actual value:
				cva::zero_floor_traits<cva::Dual<double>>().apply(dual2).value(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());

			Assert::AreEqual(
				// Expected value:
				x2 > 0 ? y2 : 0.0,
				// Actual value:
				cva::zero_floor_traits<cva::Dual<double>>().apply(dual2).deriv(),
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}
	};
}