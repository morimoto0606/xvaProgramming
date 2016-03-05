#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include <boost/bind.hpp>
#include "PathMaker.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;

	template<typename P>
	struct analytic_exposure_traits {
		typedef P payoff_type;
		typedef ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > result_type;
		static result_type apply(const double mu,
			const double sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			return result_type();
		}
	};

	template<>
	struct analytic_exposure_traits<Forward> {
		typedef Forward payoff_type;
		typedef ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > result_type;
		static result_type appry(const double mu,
			const double sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			ublas::vector<boost::function<Dual<double>(
				const Dual<double>&)> > fwdFunctions(gridNum + 1);
			for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
				const double tau = maturity
					* static_cast<double>(gridNum - gridIndex)
					/ static_cast<double>(gridNum);
				fwdFunctions(gridIndex)
					= boost::bind(forwardFunction<Dual<double> >, _1,
						Dual<double>(mu), Dual<double>(sigma),
						payoff.gearing(), payoff.strike(), tau);
			}
			return fwdFunctions;
		}
	};
	template<>
	struct analytic_exposure_traits<European> {
		typedef European payoff_type;
		typedef ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > result_type;
		static result_type appry(const double mu,
			const double sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			result_type eurFunctions(gridNum + 1);
			for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
				const double tau = maturity
					* static_cast<double>(gridNum - gridIndex)
					/ static_cast<double>(gridNum);
				eurFunctions(gridIndex)
					= boost::bind(europeanFunction<Dual<double> >, _1,
						Dual<double>(mu), Dual<double>(sigma),
						payoff.gearing(), payoff.strike(), payoff.shiftAmount(),
						tau);
			}
			return eurFunctions;
		}
	};
	template<>
	struct analytic_exposure_traits<Mountain> {
		typedef Mountain payoff_type;
		typedef ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > result_type;
		static result_type appry(const double mu,
			const double sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			result_type functions(gridNum + 1);
			for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
				const double tau = maturity
					* static_cast<double>(gridNum - gridIndex)
					/ static_cast<double>(gridNum);
				functions(gridIndex)
					= boost::bind(mountain<Dual<double> >, _1,
						Dual<double>(mu), Dual<double>(sigma),
						payoff.gearing(), payoff.strikes(),
						payoff.shiftAmount(), tau);
			}
			return functions;
		}
	};
}//namespace cva
