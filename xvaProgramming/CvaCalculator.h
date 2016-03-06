#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include "Regression.h"
#include "LsmFunctions.h"
#include <boost/bind.hpp>
#include "PathMaker.h"

namespace cva {
	namespace ublas = boost::numeric::ublas;

	template<typename Derived>
	class CvaCalculator {
	public:
		const Derived& operator()() const 
		{
			return static_cast<const Derived&>(*this);
		}
	};

	class ExplicitCalculator : public CvaCalculator<ExplicitCalculator> {
	public:
		ExplicitCalculator() {}
		template <typename T, typename P>
		T operator()(const ublas::vector<boost::function<
			T(const T&) > >& exposureFunctions,
			const Path<T>& path, const PayOff<P>& payoff) const
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
			return cvaValue * path.dt() / static_cast<double>(path.pathNum());
		}
	};

	class ImplicitCalculator : public CvaCalculator<ImplicitCalculator> {
	public:
		ImplicitCalculator() {}
		template <typename T, typename P>
		T operator()(const ublas::vector<boost::function<
			T(const T&) > >& exposureFunctions,
			const Path<T>& path, const PayOff<P>& payoff) const
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
			return cvaValue * path.dt() / static_cast<double>(path.pathNum());
		}
	};

}//namespace cva
