#include "stdafx.h"
#include "CppUnitTest.h"
#include <Payoff.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	namespace ublas = boost::numeric::ublas;
	TEST_CLASS(Payoff)
	{
	public:
		TEST_METHOD(ForwardPayoff)
		{
			ublas::vector<double> x(2);
			x(0) = 100;
			x(1) = 125.54334;
			cva::Forward payoff(2.0, 100);
			Assert::AreEqual(2.0 * (x(1) - 100), payoff(x), 1e-10);
			}
		TEST_METHOD(EuropeanPayoff)
		{
			ublas::vector<double> x(2);
			x(0) = 100;
			x(1) = 125.54334;
			cva::European payoff(2.0, 100.0, -10.0);
			Assert::AreEqual(
				// Expected value:
				2.0 * std::max(x(1) - 100.0, 0.0) - 10.0,
				// Actual value:
				payoff(x),
				// Tolerance:
				1e-5,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}

	};
}