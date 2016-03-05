#pragma once
#include "CvaCalculator.h"
namespace cva {
	template <typename C>
	struct cva_calculator_traits {
	};

	template <>
	struct cva_calculator_traits<ExplicitCalculator> {
	public:
		template <typename T>
		static const T apply()
		{
			ExplicitCalculator()
		}
	};

	template <>
	struct cva_calculator_traits<ImplicitCalculator> {
	public:
		template <typename T>
		static const T apply()
		{
			ExplicitCalculator()
		}
	};
}