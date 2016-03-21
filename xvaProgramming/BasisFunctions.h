#pragma once
#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "Path.h"

namespace ublas = boost::numeric::ublas;
namespace cva {

	template <typename T, typename U>
	struct basis_function_traits {
		typedef T value_type;
		typedef ublas::vector<U> result_type;
		typedef boost::function<result_type(const value_type&)> function_type;

		template<typename P>
		static result_type apply(const ublas::vector<const function_type>& basisFunction,
			const Path<P>& path, const size_t pathIndex, const std::size_t gridIndex)
		{
			result_type result(basisFunction.size());
			for (std::size_t i = 0; i < basisFunction.size(); ++i) {
				result(i) = basisFunction(i)(path.getPathValue(pathIndex, gridIndex));
			}
			return result;
		}
	};

	template <typename T>
	struct basis_function_traits<ublas::vector_expression<T>, T> {
		typedef ublas::vector_expression<T> value_type;
		typedef ublas::vector<T> result_type;
		typedef boost::function<result_type(const value_type&)> function_type;

		template<typename P>
		static result_type apply(const ublas::vector<const function_type>& basisFunction,
			const Path<P>& path, const std::size_t pathIndex, const std::size_t gridIndex)
		{
			result_type result(basisFunction.size());
			for (std::size_t i = 0; i < basisFunction.size(); ++i) {
				result(i) = basisFunction(i)(path.getTimewisePath(pathIndex));
			}
			return result;
		}
	};

	template <typename T, typename U>
	class BasisFunctions {
	public:
		typedef T value_type;
		typedef U result_type;
		typedef boost::function<result_type(const value_type&)> function_type;
	public:
		BasisFunctions(const ublas::vector<function_type>& basisFunctions)
			: _basisFunctions(basisFunctions) {}
		template <typename P>
		typename basis_function_traits<value_type, result_type>::result_type
			operator()(const Path<P>& path) {
			return basis_function_traits<value_type, result_type>::apply(_basisFunctions, path);
		}
	private:
		ublas::vector<function_type> _basisFunctions;
	};
}