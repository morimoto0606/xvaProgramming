#pragma once
#pragma warning(disable:4996)

double __stdcall  calcCvaLsmExposureEuropean(
	const double spot,
	const double strike,
	const double vol,
	const double maturity,
	const int gridNum,
	const long seed,
	const long pathNumForRegression,
	const long pathNumForMonte,
	const int orderOfBasisFunction,
	const bool isExplicit,
	const bool isCoeffShock);