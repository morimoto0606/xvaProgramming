#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include "Path.h"
#include "analytic_exposure_traits.h"
#include "BasisFunctions.h"
#include "Regressor.h"
#include <boost/numeric/ublas/io.hpp>

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

	template <typename T, typename R, typename S, typename P>
	class RegressionExposure : public Exposure<RegressionExposure<T, R, S, P>> {
	public:
		typedef T value_type;
		typedef T result_type;
		typedef R regression_type;
		typedef S state_type;
		typedef P payoff_type;
		typedef BasisFunctions<state_type> basis_type;
		typedef ublas::vector<value_type> coefficints_type;

		RegressionExposure(
			const Path<value_type>& path, 
			const ublas::vector<basis_type>& basisSeries,
			const Regressor<regression_type, payoff_type>& regressor)
		: _basisSeries(basisSeries),  _regressor(regressor)
		{
		}
		result_type operator()(const Path<value_type>& path, 
			const std::size_t pathIndex,
			const std::size_t gridIndx) const
		{
			BasisFunctions<state_type>::result_type basisValues
				= _basisSeries(gridIndx)(path, pathIndex, gridIndx);
			T result = ublas::inner_prod(_regressor.getCoeffs(gridIndx), basisValues);
			return result;
		}
	private:
		ublas::vector<basis_type> _basisSeries;
		Regressor<regression_type, payoff_type> _regressor;
	};
}