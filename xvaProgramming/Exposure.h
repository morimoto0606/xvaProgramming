#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include "Path.h"
#include "analytic_exposure_traits.h"
#include "BasisFunctions.h"
#include "Regression.h"

namespace ublas = boost::numeric::ublas;
namespace cva {
	template <typename Derived>
	class Exposure {
	public:
		const Derived& operator()() const
		{
			return static_cast<const Derived&>(*this);
		}
	};

	template <typename P>
	class AnalyticExposure : public Exposure<AnalyticExposure<P>> {
	public:
		typedef P payoff_type;
		AnalyticExposure(const payoff_type& payoff) : _payoff(payoff) {}
		template <typename T>
		T operator()(const Path<T>& path, const std::size_t pathIndex,
			const std::size_t gridIndx) const
		{
			return analytic_exposure_traits<payoff_type>::apply(
				path, pathIndex, gridIndx, _payoff);
		}
	private:
		payoff_type _payoff;
	};

	template <typename V, typename D, typename R, typename P>
	class RegressionExposure : public Exposure<RegressionExposure<V, D, R, P>> {
	public:
		typedef V value_type;
		typedef BasisFunctions<D, R> basis_type;
		typedef P payoff_type;
		typedef ublas::vector<value_type> coefficints_type;

		RegressionExposure(
			const Path<value_type>& path, 
			const ublas::vector<const basis_type>& basisSeries,
			const payoff_type& payoff)
		: _basisSeries(basisSeries), _payoff(payoff)
		{
			for (std::size_t gridIndex = 0; gridIndex <= path.gridNum(); ++gridIndex) {
				_coeffSeries(gridIndex) = regresssion(gridIndex, payoff, path, _basisSeries(gridIndex));
			}
		}
		value_type operator()(const Path<value_type>& path, const std::size_t pathIndex,
			const std::size_t gridIndx) const
		{
			ublas::vector<value_type> basis 
				= basis_function_traits<D, R>::apply(_basisSeries(gridIndx), path, pathIndex, gridIndx);
			T result = ublas::inner_prod(_coeffSeries(gridIndx), basis);
			return result;
		}
	private:
		payoff_type _payoff;
		ublas::vector<const basis_type> _basisSeries;
		ublas::vector<const coefficints_type> _coeffSeries;
	};
}