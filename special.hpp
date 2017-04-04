#pragma once

#include <cmath>

_STD_BEGIN

_Check_return_ double assoc_laguerre(_In_ unsigned _Degree, _In_ unsigned _Order, _In_ double _Value);
_Check_return_ inline float assoc_laguerref(_In_ unsigned _Degree, _In_ unsigned _Order, _In_ float _Value)
	{
	return (static_cast<float>(_STD assoc_laguerre(_Degree, _Order, _Value))); // narrowing
	}
_Check_return_ inline long double assoc_laguerrel(_In_ unsigned _Degree, _In_ unsigned _Order, _In_ long double _Value)
	{
	return (_STD assoc_laguerre(_Degree, _Order, static_cast<double>(_Value)));
	}

_Check_return_ double assoc_legendre(_In_ unsigned _Degree, _In_ unsigned _Order, _In_ double _Value);
_Check_return_ inline float assoc_legendref(_In_ unsigned _Degree, _In_ unsigned _Order, _In_ float _Value)
	{
	return (static_cast<float>(_STD assoc_legendre(_Degree, _Order, _Value))); // narrowing
	}
_Check_return_ inline long double assoc_legendrel(_In_ unsigned _Degree, _In_ unsigned _Order, _In_ long double _Value)
	{
	return (_STD assoc_legendre(_Degree, _Order, static_cast<double>(_Value)));
	}

_Check_return_ double legendre(_In_ unsigned _Degree, _In_ double _Value);
_Check_return_ inline float legendref(_In_ unsigned _Degree, _In_ float _Value)
	{
	return (static_cast<float>(_STD legendre(_Degree, _Value))); // narrowing
	}
_Check_return_ inline long double legendrel(_In_ unsigned _Degree, _In_ long double _Value)
	{
	return (_STD legendre(_Degree, static_cast<double>(_Value)));
	}

_Check_return_ double beta(_In_ double _Arg1, _In_ double _Arg2);
_Check_return_ inline float betaf(_In_ float _Arg1, _In_ float _Arg2)
	{
	// FIXME: narrowing:
	return (static_cast<float>(_STD beta(static_cast<double>(_Arg1), static_cast<double>(_Arg2))));
	}
_Check_return_ inline long double betal(_In_ long double _Arg1, _In_ long double _Arg2)
	{
	return (_STD beta(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_Check_return_ double hypot(_In_ double _Dx, _In_ double _Dy, _In_ double _Dz) _NOEXCEPT;
_Check_return_ inline float hypot(_In_ float const _Dx, _In_ float const _Dy, _In_ float const _Dz) _NOEXCEPT
	{
	// FIXME: narrowing:
	return (static_cast<float>(hypot(static_cast<double>(_Dx), static_cast<double>(_Dy), static_cast<double>(_Dz))));
	}
_Check_return_ inline long double hypot(_In_ long double const _Dx,
										_In_ long double const _Dy,
										_In_ long double const _Dz) _NOEXCEPT
	{
	return (hypot(static_cast<double>(_Dx), static_cast<double>(_Dy), static_cast<double>(_Dz)));
	}

_STD_END
