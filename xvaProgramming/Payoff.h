#ifndef PAYOFF_H_INCLUDED
#define PAYOFF_H_INCLUDED

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "Functions.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;
	template <typename Derived>
	class PayOff {
	public:
		const Derived& operator()() const
		{
			return static_cast<const Derived&>(*this);
		}
	};

	class Forward : public PayOff <Forward> {
	public:
		Forward(const double& a, const double& b) : _a(a), _b(b) {}
		
		template <typename T>
		typename T::value_type operator()(
			const ublas::vector_expression<T>& x) const
		{
			return _a * (*(x().end() - 1) - _b);
		}
		double gearing() const { return _a; }
		double strike() const { return _b; }
		ublas::vector<double> strikes() const {
			ublas::vector<double> x(1, _b);
			return x;
		}
		double shiftAmount() const { return 0; }
	private:
		double _a;
		double _b;
	};


	class European : public PayOff <European> {
	public:
		European(const double& a, const double& b,
			const double& c) : _a(a), _b(b), _c(c) {}
		European(const double& a, const double& b)
			: _a(a), _b(b), _c(0.0) {}
		template <typename T>
		typename T::value_type operator()(
			const ublas::vector_expression<T>& x) const
		{
			return _a * cva::zeroFloor(*(x().end() - 1) - _b) + _c;
		}
		double gearing() const { return _a; }
		double strike() const { return _b; }
		double shiftAmount() const { return _c; }
		ublas::vector<double> strikes() const {
			ublas::vector<double> x(1, _b);
			return x;
		}
	private:
		double _a;
		double _b;
		double _c;
	};

	class RiskReversal : public PayOff <RiskReversal> {
	public:
		RiskReversal(const double& a, const double& b1,
			const double& b2) 
		: _a(a), _b1(b1), _b2(b2) {}

		template <typename T>
		typename T::value_type operator()(
			const ublas::vector_expression<T>& x) const
		{
			return cva::zeroFloor(_a * *(x().end() - 1) - _b1)
				- cva::zeroFloor(_b2 -_a * *(x().end() - 1));
		}
		double gearing() const { return _a; }
		double strike1() const { return _b1; }
		double strike2() const { return _b2; }

	private:
		double _a;
		double _b1;
		double _b2;
	};

	class Mountain : public PayOff <Mountain> {
	public:
		Mountain(const double a, const double b1,
			const double b2, const double b3, const double b4)
			: _a(a), _b1(b1), _b2(b2), _b3(b3), _b4(b4), _c(0.0) {}
	
		Mountain(const double a, const double b1,
			const double b2, const double b3, const double b4,
			const double c)
		: _a(a), _b1(b1), _b2(b2), _b3(b3), _b4(b4), _c(c) {}

		template <typename T>
		typename T::value_type operator()(
			const ublas::vector_expression<T>& x) const
		{
			return cva::zeroFloor(_a * *(x().end() - 1) - _b1)
				- cva::zeroFloor(_a * *(x().end() - 1) - _b2)
				- cva::zeroFloor(_a * *(x().end() - 1) - _b3)
				+ cva::zeroFloor(_a * *(x().end() - 1) - _b4)
				+ _c;
		}

		double gearing() const { return _a; }
		double strike() const { return _b1; }
		ublas::vector<double> strikes() const {
			ublas::vector<double> x(4, 0);
			x(0) = _b1;
			x(1) = _b2;
			x(2) = _b3;
			x(3) = _b4;
			return x;
		}
		double shiftAmount() const { return _c; }
	private:
		double _a;
		double _b1;
		double _b2;
		double _b3;
		double _b4;
		double _c;
	};

	class Asian : public PayOff<Asian> {
	public:
		Asian(const double a, const double b) : _a(a), _b(b) {}
		
		template <typename T>
		typename T::value_type operator()(
			const ublas::vector_expression<T>& x) const
		{
			std::size_t gridNum = x().size() - 1;
			TimewiseAverage average{ gridNum, 1 };
			return _a * average(x) - _b;
		}

		double gearing() const { return _a; }
		double strike() const { return _b; }

	private:
		double _a;
		double _b;
	};

	class Trf : public PayOff<Trf> {
	public:
		Trf(const double a, const double b, const double c)
			: _a(a), _b(b), _c(c) {}

		template <typename T>
		typename T::value_type operator()(
			const ublas::vector_expression<T>& x) const
		{
			std::size_t gridNum = x().size() - 1;
			TimewiseAverage average{ gridNum, 1 };
			return (_a * *(x().end() - 1) - _b) * sigmoid(average(x()) - _c, 10);
		}

		double gearing() const { return _a; }
		double strike() const { return _b; }

	private:
		double _a;
		double _b;
		double _c;
	};
} // namespace cva

#endif