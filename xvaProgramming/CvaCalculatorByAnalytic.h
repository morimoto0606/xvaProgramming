#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include "Regression.h"
#include "LsmFunction.h"
#include <boost/bind.hpp>
#include "PathMaker.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;

	Dual<double> calcCvaFwdByAnalytic(
		const double x0, const double mu, const double sigma,
		const Forward& payoff, const double maturity,
		const std::size_t gridNum, const shockTypeEnum shockType)
	{
		Dual<double> cvaValue(0.0);
		Dual<double> x0Dual = shockType == undEnum
			? Dual<double>(x0, 1.0)
			: Dual<double>(x0);
		Dual<double> sigmaDual = shockType == volEnum
			? Dual<double>(sigma, 1.0)
			: Dual<double>(sigma);
		const double strike = payoff.strike();
		const double dt = maturity / static_cast<double>(gridNum);
		for (std::size_t gridIndex = 1; gridIndex <= gridNum; ++gridIndex) {
			const double t = gridIndex * dt;
			const double tau = maturity - t;
			const double gearing = std::exp(mu * tau) * payoff.gearing();
			cvaValue += forwardFunction<Dual<double> >(x0Dual,
				Dual<double>(mu), sigmaDual, gearing, strike, t);
		}
		return cvaValue;
	}

	Dual<double> calcCvaEurByAnalytic(
		const double x0, const double mu, const double sigma,
		const European& payoff, const double maturity,
		const std::size_t gridNum, const shockTypeEnum shockType)
	{
		Dual<double> cvaValue(0.0);
		Dual<double> x0Dual = shockType == undEnum
			? Dual<double>(x0, 1.0)
			: Dual<double>(x0);
		Dual<double> sigmaDual = shockType == volEnum
			? Dual<double>(sigma, 1.0)
			: Dual<double>(sigma);

		const double strike = payoff.strike();
		const double shift = payoff.shiftAmount();
		const double dt = maturity / static_cast<double>(gridNum);
		for (std::size_t gridIndex = 1; gridIndex <= gridNum; ++gridIndex) {
			const double gearing = payoff.gearing();
			const double t = gridIndex * dt;
			cvaValue += europeanFunction<Dual<double> >(x0Dual,
				Dual<double>(mu), sigmaDual, gearing, strike, shift, maturity);
		}
		return cvaValue * dt;
	}

	Dual<double> calcCvaMountainByAnalytic(
		const double x0, const double mu, const double sigma,
		const Mountain& payoff, const double maturity,
		const std::size_t gridNum, const shockTypeEnum shockType)
	{
		Dual<double> cvaValue(0.0);
		Dual<double> x0Dual = shockType == undEnum
			? Dual<double>(x0, 1.0)
			: Dual<double>(x0);
		Dual<double> sigmaDual = shockType == volEnum
			? Dual<double>(sigma, 1.0)
			: Dual<double>(sigma);

		const double dt = maturity / static_cast<double>(gridNum);
		for (std::size_t gridIndex = 1; gridIndex <= gridNum; ++gridIndex) {
			const double gearing = payoff.gearing();
			const double t = gridIndex * dt;
			cvaValue += mountain<Dual<double> >(x0Dual,
				Dual<double>(mu), sigmaDual, gearing, payoff.strikes(),
				payoff.shiftAmount(), maturity);
		}
		return cvaValue * dt;
	}
}//namespace cva
