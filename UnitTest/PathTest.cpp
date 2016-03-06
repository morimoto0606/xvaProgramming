#include "stdafx.h"
#include "CppUnitTest.h"
#include <Path.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace cva;
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
			std::size_t seed = 100;
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
		TEST_METHOD(PathStatistics)
		{
			const Dual<double> x0(100.0, 1.0);
			const Dual<double> mu(0.2);
			const Dual<double> sigma(0.5);
			const std::size_t pathNum = 100000;
			const std::size_t gridNum = 1;
			const double dt = 1;
			const std::size_t seed = 1;
			cva::Path<Dual<double>> path(
				x0, mu, sigma, pathNum,
				gridNum, dt, seed);

			double logMean = 0.0;
			double logMeanSquare = 0.0;
			for (std::size_t i = 0; i < pathNum; ++i) {
				logMean += std::log(
					path.getPathValue(
						i, gridNum).value() / 100.0);
				logMeanSquare += std::pow(
					std::log(path.getPathValue(
						i, gridNum).value() / 100.0), 2.0);
			}
			logMean /= pathNum;
			const double var
				= logMeanSquare / pathNum - std::pow(logMean, 2.0);

			Assert::AreEqual(0.0758121, logMean, 1e-5,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
			Assert::AreEqual(0.499263, std::sqrt(var), 1e-5,
				// Message:
				L"Basic test failed",
				// Line number - used if there is no PDB file:
				LINE_INFO());
		}
	};
}