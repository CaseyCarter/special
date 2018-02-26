#pragma once

#include <cmath>
#include <type_traits>

// For AppVeyor
#ifndef _NODISCARD
#define _NODISCARD
#endif

enum { PM_TEST_PASS = 0, PM_TEST_FAIL = 1 };

_STD_BEGIN
#if _HAS_CXX17
_NODISCARD double assoc_laguerre(unsigned _Degree, unsigned _Order, double _Value);
_NODISCARD float assoc_laguerref(unsigned _Degree, unsigned _Order, float _Value);
_NODISCARD inline long double assoc_laguerrel(const unsigned _Degree, const unsigned _Order,
	const long double _Value)
	{
	return (_STD assoc_laguerre(_Degree, _Order, static_cast<double>(_Value)));
	}

_NODISCARD double assoc_legendre(unsigned _Degree, unsigned _Order, double _Value);
_NODISCARD float assoc_legendref(unsigned _Degree, unsigned _Order, float _Value);
_NODISCARD inline long double assoc_legendrel(const unsigned _Degree, const unsigned _Order,
	const long double _Value)
	{
	return (_STD assoc_legendre(_Degree, _Order, static_cast<double>(_Value)));
	}

_NODISCARD double beta(double _Arg1, double _Arg2);
_NODISCARD float betaf(float _Arg1, float _Arg2);
_NODISCARD inline long double betal(const long double _Arg1, const long double _Arg2)
	{
	return (_STD beta(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double comp_ellint_1(double);
_NODISCARD float comp_ellint_1f(float);
_NODISCARD inline long double comp_ellint_1l(const long double _Arg)
	{
	return (_STD comp_ellint_1(static_cast<double>(_Arg)));
	}

_NODISCARD double comp_ellint_2(double);
_NODISCARD float comp_ellint_2f(float);
_NODISCARD inline long double comp_ellint_2l(const long double _Arg)
	{
	return (_STD comp_ellint_2(static_cast<double>(_Arg)));
	}

_NODISCARD double comp_ellint_3(double, double);
_NODISCARD float comp_ellint_3f(float, float);
_NODISCARD inline long double comp_ellint_3l(const long double _Arg1, const long double _Arg2)
	{
	return (_STD comp_ellint_3(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_bessel_i(double, double);
_NODISCARD float cyl_bessel_if(float, float);
_NODISCARD inline long double cyl_bessel_il(const long double _Arg1, const long double _Arg2)
	{
	return (_STD cyl_bessel_i(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_bessel_j(double, double);
_NODISCARD float cyl_bessel_jf(float, float);
_NODISCARD inline long double cyl_bessel_jl(const long double _Arg1, const long double _Arg2)
	{
	return (_STD cyl_bessel_j(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_bessel_k(double, double);
_NODISCARD float cyl_bessel_kf(float, float);
_NODISCARD inline long double cyl_bessel_kl(const long double _Arg1, const long double _Arg2)
	{
	return (_STD cyl_bessel_k(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_neumann(double, double);
_NODISCARD float cyl_neumannf(float, float);
_NODISCARD inline long double cyl_neumannl(const long double _Arg1, const long double _Arg2)
	{
	return (_STD cyl_neumann(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double ellint_1(double, double);
_NODISCARD float ellint_1f(float, float);
_NODISCARD inline long double ellint_1l(const long double _Arg1, const long double _Arg2)
	{
	return (_STD ellint_1(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double ellint_2(double, double);
_NODISCARD float ellint_2f(float, float);
_NODISCARD inline long double ellint_2l(const long double _Arg1, const long double _Arg2)
	{
	return (_STD ellint_2(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double ellint_3(double, double, double);
_NODISCARD float ellint_3f(float, float, float);
_NODISCARD inline long double ellint_3l(const long double _Arg1, const long double _Arg2,
	const long double _Arg3)
	{
	return (_STD ellint_3(static_cast<double>(_Arg1), static_cast<double>(_Arg2),
		static_cast<double>(_Arg3)));
	}

_NODISCARD double expint(double);
_NODISCARD float expintf(float);
_NODISCARD inline long double expintl(const long double _Arg)
	{
	return (_STD expint(static_cast<double>(_Arg)));
	}

_NODISCARD double hermite(unsigned, double);
_NODISCARD float hermitef(unsigned, float);
_NODISCARD inline long double hermitel(const unsigned _Arg1, const long double _Arg2)
	{
	return (_STD hermite(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double laguerre(unsigned, double);
_NODISCARD float laguerref(unsigned, float);
_NODISCARD inline long double laguerrel(const unsigned _Arg1, const long double _Arg2)
	{
	return (_STD laguerre(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double legendre(unsigned _Degree, double _Value);
_NODISCARD float legendref(unsigned _Degree, float _Value);
_NODISCARD inline long double legendrel(const unsigned _Degree, const long double _Value)
	{
	return (_STD legendre(_Degree, static_cast<double>(_Value)));
	}

_NODISCARD double riemann_zeta(double);
_NODISCARD float riemann_zetaf(float);
_NODISCARD inline long double riemann_zetal(const long double _Arg)
	{
	return (_STD riemann_zeta(static_cast<double>(_Arg)));
	}

_NODISCARD double sph_bessel(unsigned, double);
_NODISCARD float sph_besself(unsigned, float);
_NODISCARD inline long double sph_bessell(const unsigned _Arg1, const long double _Arg2)
	{
	return (_STD sph_bessel(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double sph_legendre(unsigned, unsigned, double);
_NODISCARD float sph_legendref(unsigned, unsigned, float);
_NODISCARD inline long double sph_legendrel(const unsigned _Arg1, const unsigned _Arg2,
	const long double _Theta)
	{
	return (_STD sph_legendre(_Arg1, _Arg2, static_cast<double>(_Theta)));
	}

_NODISCARD double sph_neumann(unsigned, double);
_NODISCARD float sph_neumannf(unsigned, float);
_NODISCARD inline long double sph_neumannl(const unsigned _Arg1, const long double _Arg2)
	{
	return (_STD sph_neumann(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double hypot(double _Dx, double _Dy, double _Dz);
_NODISCARD float hypot(float _Dx, float _Dy, float _Dz);
_NODISCARD inline long double hypot(const long double _Dx, const long double _Dy,
	const long double _Dz)
	{
	return (_STD hypot(static_cast<double>(_Dx), static_cast<double>(_Dy), static_cast<double>(_Dz)));
	}

template<class _Ty1,
	class _Ty2,
	class _Ty3,
	enable_if_t<is_arithmetic_v<_Ty1> && is_arithmetic_v<_Ty2>
		&& is_arithmetic_v<_Ty3>, int> = 0>
	_NODISCARD inline auto hypot(const _Ty1 _Dx, const _Ty2 _Dy, const _Ty3 _Dz)
	{	// bring mixed types to a common type
	using _Common = _Common_float_type_t<_Ty1, _Common_float_type_t<_Ty2, _Ty3>>;
	return (_STD hypot(static_cast<_Common>(_Dx), static_cast<_Common>(_Dy),
		static_cast<_Common>(_Dz)));
	}
#endif /* _HAS_CXX17 */
_STD_END
