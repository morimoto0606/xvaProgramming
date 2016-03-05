#include "stdafx.h"
#include "CppUnitTest.h"
#include <Path.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	namespace ublas = boost::numeric::ublas;
	TEST_CLASS(Path)
	{
	public:
		TEST_METHOD(PathSize)
		{
			ublas::vector<double> x(2);
			double x0 = 100;
			double mu = 0.1;
			double sigma = 0.3;
			std::size_t pathNum = 1;
			std::size_t gridNum = 10;
			double dt = 0.1;
			double seed = 100;
			cva::Path<double> 
				path(x0, mu, sigma, pathNum, gridNum, dt, seed);

			Assert::AreEqual(
				// Expected value:
				path.pathNum(),
				// Actual value:
				pathNum,
				// Tolerance:
				1e-10,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}
	};
}