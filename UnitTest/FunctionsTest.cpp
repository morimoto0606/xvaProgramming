#include "stdafx.h"
#include "CppUnitTest.h"
#include <Functions.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace cva;
namespace UnitTest
{
	namespace ublas = boost::numeric::ublas;
	TEST_CLASS(Functions)
	{
	public:
		TEST_METHOD(PathwiseMonomial)
		{
			ublas::vector<double> x(5);
			x(0) = 1212;
			x(1) = 34378;
			x(2) = 0.3433;
			x(3) = -34343.3243;
			x(4) = 3432379872.3234;
			cva::PathwiseMonomial foo(3, 1);
			double d = foo(x);
			Assert::AreEqual(
				x(3), d, 1e-10);
		}
		TEST_METHOD(PathwiseSum)
		{
			ublas::vector<double> x(5);
			x(0) = 1212.0;
			x(1) = 34378.0;
			x(2) = 0.3433;
			x(3) = -34343.3243;
			x(4) = 3432379872.3234;
			cva::PathwiseSum foo(3, 1);
			double d = foo(x);
			double expect = x(0) + x(1) + x(2) + x(3);
			Assert::AreEqual(
				// Expected value:
				expect,
				// Actual value:
				d,
				// Tolerance:
				1e-5,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}

	};
}