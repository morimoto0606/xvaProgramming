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

	template<typename T, typename U>
	class CvaCalculator {
	public:
	private:

	};
	template <typename T>
	T calcCvaUsingExplicitExposure(
		const ublas::vector<boost::function<
		T(const T&) > >& exposureFunctions,
		const Path<T>& path, const double dt)
	{
		T cvaValue(0.0);
		for (std::size_t pathIndex = 0; 
		pathIndex < path.pathNum(); ++pathIndex) {
			T pathwiseValue(0.0);
			for (std::size_t gridIndex = 1; 
			gridIndex <= path.gridNum(); ++gridIndex) {
				pathwiseValue += cva::zeroFloor(
					exposureFunctions(gridIndex)(
						path.getPathValue(pathIndex, gridIndex)));
			}
			cvaValue += pathwiseValue;
		}
		return cvaValue * dt / static_cast<double>(path.pathNum());
	}

	template <typename T, typename U>
	T calcCvaUsingImplicitExposure(
		const ublas::vector<boost::function<
		T (const T&) > >& exposureFunctions,
		const Path<T>& path, const PayOff<U>& payoff,
		const double dt)
	{
		T cvaValue(0.0);
		for (std::size_t pathIndex = 0;
		pathIndex < path.pathNum(); ++pathIndex) {
			T pathwiseValue(0.0);
			for (std::size_t gridIndex = 1;
			gridIndex <= path.gridNum(); ++gridIndex) {
				if ((exposureFunctions(gridIndex)(
					path.getPathValue(pathIndex, gridIndex))).value() > 0.0) {
					pathwiseValue += 
						payoff()(path.getTimewisePath(pathIndex));
				}
			}
			cvaValue += pathwiseValue;
		}
		return cvaValue * dt / static_cast<double>(path.pathNum());
	}

}//namespace cva
