#pragma once

#include <cmath>

_STD_BEGIN

double assoc_laguerre(unsigned _Degree, unsigned _Order, double _Value);
inline float assoc_laguerref(unsigned _Degree, unsigned _Order, float _Value)
	{
	return (static_cast<float>(_STD assoc_laguerre(_Degree, _Order, _Value))); // narrowing
	}
inline long double assoc_laguerrel(unsigned _Degree, unsigned _Order, long double _Value)
	{
	return (_STD assoc_laguerre(_Degree, _Order, static_cast<double>(_Value)));
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

double legendre(unsigned _Degree, double _Value);
inline float legendref(unsigned _Degree, float _Value)
	{
	return (static_cast<float>(_STD legendre(_Degree, _Value))); // narrowing
	}
inline long double legendrel(unsigned _Degree, long double _Value)
	{
	return (_STD legendre(_Degree, static_cast<double>(_Value)));
	}

double beta(double _Arg1, double _Arg2);
inline float betaf(float _Arg1, float _Arg2)
	{
	return (static_cast<float>(_STD beta(static_cast<double>(_Arg1), static_cast<double>(_Arg2)))); // narrowing
	}
inline long double betal(long double _Arg1, long double _Arg2)
	{
	return (_STD beta(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_STD_END
