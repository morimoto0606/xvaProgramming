#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include "Path.h"
#include "BasisFunctions.h"
#include "Regression.h"

namespace ublas = boost::numeric::ublas;

namespace {
	template<typename T, typename P>
	T getPv(const cva::Path<T>& path, const cva::PayOff<P>& payoff) {
		T pv(0.0);
		for (std::size_t i = 0; i < path.pathNum(); ++i)
		{
			pv += payoff()(path.getTimewisePath(i));
		}
		return pv / static_cast<double>(path.pathNum());
	}

}
namespace cva {
	template <typename T, typename P>
	class Regressor {
	public:
		typedef T value_type;
		typedef BasisFunctions<T> basis_type;
		typedef P payoff_type;
		typedef ublas::vector<value_type> coefficints_type;

	public:
		Regressor(
			const Path<value_type>& path,
			const ublas::vector<basis_type>& basisSeries,
			const payoff_type& payoff)
			: _basisSeries(basisSeries), _payoff(payoff)
		{
			_coeffSeries = ublas::vector<coefficints_type>(path.gridNum());
			//Regression is impossible for gridIndex = 0, EE(0) = pv
			_coeffSeries(0) = ublas::zero_vector<value_type>(_basisSeries(0).size());
			_coeffSeries(0)(0) = getPv(path, payoff);
			//Regression is possible for gridIndex > 0
			for (std::size_t gridIndex = 1; gridIndex < path.gridNum(); ++gridIndex) {
				_coeffSeries(gridIndex) = regresssion(gridIndex, payoff, path, _basisSeries(gridIndex));
			}
		}

		coefficints_type getCoeffs(const std::size_t gridIndex) const
		{
			return _coeffSeries(gridIndex);
		}
	private:
		payoff_type _payoff;
		ublas::vector<basis_type> _basisSeries;
		ublas::vector<coefficints_type> _coeffSeries;
	};
}