#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include <boost/bind.hpp>

namespace cva {
	namespace ublas = boost::numeric::ublas;

	template<typename P>
	struct analytic_exposure_traits {
	public:
		typedef P payoff_type;
	};

	template<>
	struct analytic_exposure_traits<Forward> {
		typedef Forward payoff_type;

		template <typename T>
		static const T apply(const Path<T>& path, const std::size_t pathIndex,
			const std::size_t gridIndex, const payoff_type& payoff)
		{
			const std::size_t gridNum = path.gridNum();
			const double maturity = path.maturity();
			const double tau = maturity
				* static_cast<double>(gridNum - gridIndex)
				/ static_cast<double>(gridNum);
			return forwardFunction<T>(
				path.getPathValue(pathIndex, gridIndex),
				path.mu(),
				path.sigma(),
				payoff.gearing(),
				payoff.strike(),
				tau);
		}
	};

	template<>
	struct analytic_exposure_traits<European> {
		typedef European payoff_type;

		template <typename T>
		static const T apply(const Path<T>& path, const std::size_t pathIndex,
			const std::size_t gridIndex, const payoff_type& payoff)
		{
			const std::size_t gridNum = path.gridNum();
			const double maturity = path.maturity();
			const double tau = maturity
				* static_cast<double>(gridNum - gridIndex)
				/ static_cast<double>(gridNum);
			return europeanFunction<T>(
				path.getPathValue(pathIndex, gridIndex),
				path.mu(),
				path.sigma(),
				payoff.gearing(),
				payoff.strike(),
				payoff.shiftAmount(),
				tau);
		}
	};

	template<>
	struct analytic_exposure_traits<Mountain> {
		typedef Mountain payoff_type;

		template <typename T>
		static const T apply(const Path<T>& path, const std::size_t pathIndex,
			const std::size_t gridIndex, const payoff_type& payoff)
		{
			const std::size_t gridNum = path.gridNum();
			const double maturity = path.maturity();
			const double tau = maturity
				* static_cast<double>(gridNum - gridIndex)
				/ static_cast<double>(gridNum);
			return mountain<T>(path.getPathValue(pathIndex, gridIndex),
				path.mu(), path.sigma(),
				payoff.gearing(), payoff.strikes(),
				payoff.shiftAmount(), tau);
		}
	};

	template<>
	struct analytic_exposure_traits<Asian> {
		typedef Asian payoff_type;

		template <typename T>
		static const T apply(const Path<T>& path, const std::size_t pathIndex,
			const std::size_t gridIndex, const payoff_type& payoff)
		{
			const std::size_t gridNum = path.gridNum();
			const double maturity = path.maturity();
			ublas::vector<double> taus(gridNum - gridIndex);
			for (std::size_t i = 0; i < gridNum -gridIndex; ++i) {
				taus(i) = path.dt() * static_cast<double>(i + 1);
			}
			return asian<T>(path.getTimewisePath(pathIndex),
				taus, path.mu(),
				payoff.gearing(), payoff.strike(),
				gridIndex, gridNum);
		}
	};
}//namespace cva
