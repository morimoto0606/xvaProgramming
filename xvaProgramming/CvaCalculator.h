#pragma once
#include <iostream>
#include "Dual.h"
#include "AnalyticOptionFunctions.h"
#include "Path.h"
#include "PayOff.h"
#include "Regression.h"
#include "Exposure.h"
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
		template <typename T, typename P, typename E>
		T operator()(const Exposure<E>& exposure,
			const Path<T>& path, const PayOff<P>& payoff) const
		{
			T cvaValue(0.0);
			for (std::size_t pathIndex = 0;
				pathIndex < path.pathNum(); ++pathIndex) {
				T pathwiseValue(0.0);
				for (std::size_t gridIndex = 1;
					gridIndex <= path.gridNum(); ++gridIndex) {
					pathwiseValue += cva::zeroFloor(
						exposrue()(path, pathIndex, gridIndex));
				}
				cvaValue += pathwiseValue;
			}
			return cvaValue * path.dt() / static_cast<double>(path.pathNum());
		}
	};

	class ImplicitCalculator : public CvaCalculator<ImplicitCalculator> {
	public:
		ImplicitCalculator() {}
		template <typename T, typename P, typename E>
		T operator()(const Exposure<E>& exposure,
			const Path<T>& path, const PayOff<P>& payoff) const
		{
			T cvaValue(0.0);
			for (std::size_t pathIndex = 0;
			pathIndex < path.pathNum(); ++pathIndex) {
				T pathwiseValue(0.0);
				for (std::size_t gridIndex = 1;
				gridIndex <= path.gridNum(); ++gridIndex) {
					if (exposure()(path, pathIndex, gridIndex) > 0.0) {
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
