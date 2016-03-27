#pragma once
#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "Path.h"

namespace ublas = boost::numeric::ublas;
namespace cva {

	template <typename T>
	struct basis_function_traits {
		typedef T value_type;
		typedef ublas::vector<T> result_type;
		typedef boost::function<value_type(const value_type&)> function_type;

		template<typename P>
		static result_type apply(const ublas::vector<function_type>& basisFunction,
			const Path<P>& path, const size_t pathIndex, const std::size_t gridIndex)
		{
			result_type result(basisFunction.size());
			for (std::size_t i = 0; i < basisFunction.size(); ++i) {
				result(i) = (basisFunction(i))(path.getPathValue(pathIndex, gridIndex));
			}
			return result;
		}
	};

	template <typename T>
	struct basis_function_traits<ublas::vector<T>> {
		typedef ublas::vector<T> value_type;
		typedef ublas::vector<T> result_type;
		typedef boost::function<result_type(const value_type&)> function_type;

		template<typename P>
		static result_type apply(const ublas::vector<function_type>& basisFunction,
			const Path<P>& path, const std::size_t pathIndex, const std::size_t gridIndex)
		{
			result_type result(basisFunction.size());
			for (std::size_t i = 0; i < basisFunction.size(); ++i) {
				result(i) = (basisFunction(i))(path.getTimewisePath(pathIndex));
			}
			return result;
		}
	};

	template <typename T>
	class BasisFunctions {
	public:
		typedef T value_type;
		typedef typename basis_function_traits<value_type>::result_type result_type;
		typedef boost::function<value_type(const value_type&)> function_type;

		BasisFunctions() {}
		explicit BasisFunctions(const ublas::vector<function_type>& basisFunctions)
			: _basisFunctions(basisFunctions) {}
		template <typename P>
		result_type operator()(const Path<P>& path, const std::size_t pathIndex,
				const std::size_t gridIndex) const {
			return basis_function_traits<value_type>::apply(_basisFunctions, path,
				pathIndex, gridIndex);
		}
		std::size_t size() const {
			return _basisFunctions.size();
		}
	private:
		ublas::vector<function_type> _basisFunctions;
	};
}