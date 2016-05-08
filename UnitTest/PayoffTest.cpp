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
			Assert::AreEqual(2.0 * std::max(x(1) - 100.0, 0.0) - 10.0, payoff(x), 1e-5);
		}

		TEST_METHOD(AsianPayoff)
		{
			ublas::vector<double> x(6);
			x(0) = 100;
			x(1) = 125.54334;
			x(2) = -11232.1212;
			x(3) = -0.4343432424;
			x(4) = 433.43243;
			x(5) = 455.4546;
			cva::Asian payoff(1.0, 0.0);
			const double expected = std::accumulate(x.begin(), x.begin() + 6, 0.0) / 6.0;
			const double actual = payoff(x);
			Assert::AreEqual(expected, actual, 1e-5);
		}

	};
}