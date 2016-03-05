#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include <boost/function.hpp>

namespace cva {
	namespace ublas = boost::numeric::ublas;
	template <typename T, typename C>
	class LsmFunction {
	public:
		typedef T value_type;
		typedef C coeffcient_type;
		typedef std::size_t size_type;
		typedef value_type result_type;
		typedef boost::function<result_type(const value_type&)> function_type;

		LsmFunction() {}
		LsmFunction(const ublas::vector<coeffcient_type>& coeffs,
			const ublas::vector <boost::function<result_type (const value_type&)> >&
			basisFunctions) : _coefficients(coeffs), _basisFunctions(basisFunctions) {}
		
		result_type operator()(const value_type& x) const
		{
			ublas::vector<result_type> value(_coefficients.size());
			for (std::size_t i = 0; i < _coefficients.size(); ++i) {
				value(i) = (_basisFunctions(i))(x);
			}
			return ublas::inner_prod(value, _coefficients);
		}
	private:
		ublas::vector<coeffcient_type> _coefficients;
		ublas::vector <function_type>_basisFunctions;
	};
}//namespace cva