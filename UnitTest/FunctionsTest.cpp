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
		TEST_METHOD(Power)
		{
			const double x = -0.23232;
			const double y = -0.432;
			const std::size_t order = 5;
			const Dual<double> dual(x, y);
			const Dual<double> actual = power(dual, order);
			const double expectV = -std::pow(abs(x), order);
			const double expectD = -static_cast<double>(order) * std::pow(abs(x), order - 1) * abs(y);
			Assert::AreEqual(expectV, actual.value(), 1e-5);
			Assert::AreEqual(expectD, actual.deriv(), 1e-5);
		}
		TEST_METHOD(TimewiseAverage)
		{
			ublas::vector<double> x(5);
			x(0) = 2.0;
			x(1) = -8.0;
			x(2) = 0.3433;
			x(3) = -3.3243;
			x(4) = -2.3234;
			cva::TimewiseAverage foo(4, 1);
			//double
			{
				double d = foo(x);
				double expect = (x(0) + x(1) + x(2) + x(3)) / 4.0;
				Assert::AreEqual(expect, d, 1e-5);
			}
			ublas::vector<double> y(5);
			y(0) = -3.0;
			y(1) = -4.0;
			y(2) = 0.069354385784;
			y(3) = -1.3243;
			y(4) = 3.35482;
			
			//Dual<double> 
			{
				ublas::vector <cva::Dual<double>> z(5);
				for (std::size_t i = 0; i < 5; ++i) {
					z(i) = cva::Dual<double>(x(i), y(i));
				}
				Dual<double> d = foo(z);
				double expectV = (x(0) + x(1) + x(2) + x(3)) / 4.0;
				double expectD= (y(0) + y(1) + y(2) + y(3)) / 4.0;
				Assert::AreEqual(expectV, d.value(), 1e-2);
				Assert::AreEqual(expectD, d.deriv(), 1e-2);
			}
			cva::TimewiseAverage buz(5, 2);
			//Dual<double> 
			{
				ublas::vector <cva::Dual<double>> z(5);
				for (std::size_t i = 0; i < 5; ++i) {
					z(i) = cva::Dual<double>(x(i), y(i));
				}
				Dual<double> d = buz(z);
				double aveV = (x(0) + x(1) + x(2) + x(3) + x(4)) / 5.0;
				double aveD = (y(0) + y(1) + y(2) + y(3) + y(4)) / 5.0;
				double expectV = std::pow(aveV, 2);
				double expectD = 2.0 * aveV * aveD;

				Assert::AreEqual(expectV, d.value(), 1e-2);
				Assert::AreEqual(expectD, d.deriv(), 1e-2);
			}
		}
	};
}