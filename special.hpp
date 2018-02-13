#pragma once

#include <cmath>

_STD_BEGIN

_NODISCARD double assoc_laguerre(unsigned _Degree, unsigned _Order, double _Value);
_NODISCARD float assoc_laguerref(unsigned _Degree, unsigned _Order, float _Value);
_NODISCARD inline long double assoc_laguerrel(unsigned _Degree, unsigned _Order, long double _Value)
	{
	return (_STD assoc_laguerre(_Degree, _Order, static_cast<double>(_Value)));
	}

_NODISCARD double assoc_legendre(unsigned _Degree, unsigned _Order, double _Value);
_NODISCARD float assoc_legendref(unsigned _Degree, unsigned _Order, float _Value);
_NODISCARD inline long double assoc_legendrel(unsigned _Degree, unsigned _Order, long double _Value)
	{
	return (_STD assoc_legendre(_Degree, _Order, static_cast<double>(_Value)));
	}

_NODISCARD double legendre(unsigned _Degree, double _Value);
_NODISCARD float legendref(unsigned _Degree, float _Value);
_NODISCARD inline long double legendrel(unsigned _Degree, long double _Value)
	{
	return (_STD legendre(_Degree, static_cast<double>(_Value)));
	}

_NODISCARD double beta(double _Arg1, double _Arg2);
_NODISCARD float betaf(float _Arg1, float _Arg2);
_NODISCARD inline long double betal(long double _Arg1, long double _Arg2)
	{
	return (_STD beta(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double hypot(double _Dx, double _Dy, double _Dz) _NOEXCEPT;
_NODISCARD float hypot(float const _Dx, float const _Dy, float const _Dz) _NOEXCEPT;
_NODISCARD inline long double hypot(long double const _Dx, long double const _Dy,
	long double const _Dz) _NOEXCEPT
	{
	return (hypot(static_cast<double>(_Dx), static_cast<double>(_Dy), static_cast<double>(_Dz)));
	}

_STD_END
