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
	public:
		typedef P payoff_type;
	public:
		template <typename T>
		static ublas::vector<boost::function<T (const T&)>> apply(
			const T& mu, const T& sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			return ublas::vector<boost::function<T(const T&)>>();
		}
	};

	template<>
	struct analytic_exposure_traits<Forward> {
		typedef Forward payoff_type;
		typedef ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > result_type;
		template <typename T>
		static ublas::vector<boost::function<T(const T&)>> apply(
			const T& mu, const T& sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			ublas::vector<boost::function<T (const T&)>> fwdFunctions(gridNum + 1);
			for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
				const double tau = maturity
					* static_cast<double>(gridNum - gridIndex)
					/ static_cast<double>(gridNum);
				fwdFunctions(gridIndex)
					= boost::bind(forwardFunction<T>, _1,
						mu, sigma, payoff.gearing(), payoff.strike(), tau);
			}
			return fwdFunctions;
		}
	};
	template<>
	struct analytic_exposure_traits<European> {
		typedef European payoff_type;
		typedef ublas::vector<boost::function<Dual<double>(
			const Dual<double>&)> > result_type;
		
		template <typename T>
		static ublas::vector<boost::function<T(const T&)>> apply(
			const T& mu, const T& sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			ublas::vector<boost::function<T(const T&)>> eurFunctions(gridNum + 1);
			for (std::size_t gridIndex = 0; gridIndex <= gridNum; ++gridIndex) {
				const double tau = maturity
					* static_cast<double>(gridNum - gridIndex)
					/ static_cast<double>(gridNum);
				eurFunctions(gridIndex)
					= boost::bind(europeanFunction<T>, _1,
						mu, sigma, payoff.gearing(), payoff.strike(), payoff.shiftAmount(),
						tau);
			}
			return eurFunctions;
		}
	};
	template<>
	struct analytic_exposure_traits<Mountain> {
		typedef Mountain payoff_type;
		template <typename T>
		static ublas::vector<boost::function<T(const T&)>> apply(
			const T& mu, const T& sigma, const std::size_t gridNum,
			const double maturity, const payoff_type& payoff)
		{
			ublas::vector<boost::function<T(const T&)>> functions(gridNum + 1);
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
