#pragma once

#include <cmath>

_STD_BEGIN

double legendre(unsigned _Degree, double _Value);
inline float legendref(unsigned _Degree, float _Value)
	{
	return (static_cast<float>(_STD legendre(_Degree, _Value))); // narrowing
	}
inline long double legendrel(unsigned _Degree, long double _Value)
	{
	return (_STD legendre(_Degree, static_cast<double>(_Value)));
	}

double assoc_legendre(unsigned _Degree, unsigned _Order, double _Value);
inline float assoc_legendref(unsigned _Degree, unsigned _Order, float _Value)
	{
	return (static_cast<float>(_STD assoc_legendre(_Degree, _Order, _Value))); // narrowing
	}
inline long double assoc_legendrel(unsigned _Degree, unsigned _Order, long double _Value)
	{
	return (_STD assoc_legendre(_Degree, _Order, static_cast<double>(_Value)));
	}

double assoc_laguerre(unsigned _Degree, unsigned _Order, double _Value);
inline float assoc_laguerref(unsigned _Degree, unsigned _Order, float _Value)
	{
	return (static_cast<float>(_STD assoc_laguerre(_Degree, _Order, _Value))); // narrowing
	}
inline long double assoc_laguerrel(unsigned _Degree, unsigned _Order, long double _Value)
	{
	return (_STD assoc_laguerre(_Degree, _Order, static_cast<double>(_Value)));
	}

_STD_END
