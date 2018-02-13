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

_NODISCARD double beta(double _Arg1, double _Arg2);
_NODISCARD float betaf(float _Arg1, float _Arg2);
_NODISCARD inline long double betal(long double _Arg1, long double _Arg2)
	{
	return (_STD beta(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double comp_ellint_1(double);
_NODISCARD float comp_ellint_1f(float);
_NODISCARD inline long double comp_ellint_1l(long double _Arg)
	{
	return (_STD comp_ellint_1(static_cast<double>(_Arg)));
	}

_NODISCARD double comp_ellint_2(double);
_NODISCARD float comp_ellint_2f(float);
_NODISCARD inline long double comp_ellint_2l(long double _Arg)
	{
	return (_STD comp_ellint_2(static_cast<double>(_Arg)));
	}

_NODISCARD double comp_ellint_3(double, double);
_NODISCARD float comp_ellint_3f(float, float);
_NODISCARD inline long double comp_ellint_3l(long double _Arg1, long double _Arg2)
	{
	return (_STD comp_ellint_3(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_bessel_i(double, double);
_NODISCARD float cyl_bessel_if(float, float);
_NODISCARD inline long double cyl_bessel_il(long double _Arg1, long double _Arg2)
	{
	return (_STD cyl_bessel_i(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_bessel_j(double, double);
_NODISCARD float cyl_bessel_jf(float, float);
_NODISCARD inline long double cyl_bessel_jl(long double _Arg1, long double _Arg2)
	{
	return (_STD cyl_bessel_j(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_bessel_k(double, double);
_NODISCARD float cyl_bessel_kf(float, float);
_NODISCARD inline long double cyl_bessel_kl(long double _Arg1, long double _Arg2)
	{
	return (_STD cyl_bessel_k(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double cyl_neumann(double, double);
_NODISCARD float cyl_neumannf(float, float);
_NODISCARD inline long double cyl_neumannl(long double _Arg1, long double _Arg2)
	{
	return (_STD cyl_neumann(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double ellint_1(double k, double phi);
_NODISCARD float ellint_1f(float k, float phi);
_NODISCARD inline long double ellint_1l(long double _Arg1, long double _Arg2)
	{
	return (_STD ellint_1(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double ellint_2(double k, double phi);
_NODISCARD float ellint_2f(float k, float phi);
_NODISCARD inline long double ellint_2l(long double _Arg1, long double _Arg2)
	{
	return (_STD ellint_2(static_cast<double>(_Arg1), static_cast<double>(_Arg2)));
	}

_NODISCARD double ellint_3(double k, double, double phi);
_NODISCARD float ellint_3f(float k, float, float phi);
_NODISCARD inline long double ellint_3l(long double _Arg1, long double _Arg2, long double _Arg3)
	{
	return (_STD ellint_3(static_cast<double>(_Arg1), static_cast<double>(_Arg2),
		static_cast<double>(_Arg3)));
	}

_NODISCARD double expint(double);
_NODISCARD float expintf(float);
_NODISCARD inline long double expintl(long double _Arg)
	{
	return (_STD expint(static_cast<double>(_Arg)));
	}

_NODISCARD double hermite(unsigned n, double);
_NODISCARD float hermitef(unsigned n, float);
_NODISCARD inline long double hermitel(unsigned _Arg1, long double _Arg2)
	{
	return (_STD hermite(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double laguerre(unsigned n, double);
_NODISCARD float laguerref(unsigned n, float);
_NODISCARD inline long double laguerrel(unsigned _Arg1, long double _Arg2)
	{
	return (_STD laguerre(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double legendre(unsigned _Degree, double _Value);
_NODISCARD float legendref(unsigned _Degree, float _Value);
_NODISCARD inline long double legendrel(unsigned _Degree, long double _Value)
	{
	return (_STD legendre(_Degree, static_cast<double>(_Value)));
	}

_NODISCARD double riemann_zeta(double);
_NODISCARD float riemann_zetaf(float);
_NODISCARD inline long double riemann_zetal(long double _Arg)
	{
	return (_STD riemann_zeta(_Arg));
	}

_NODISCARD double sph_bessel(unsigned, double);
_NODISCARD float sph_besself(unsigned, float);
_NODISCARD inline long double sph_bessell(unsigned _Arg1, long double _Arg2)
	{
	return (_STD sph_bessel(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double sph_legendre(unsigned, unsigned, double);
_NODISCARD float sph_legendref(unsigned, unsigned, float);
_NODISCARD inline long double sph_legendrel(unsigned _Arg1, unsigned _Arg2, long double _Theta)
	{
	return (_STD sph_legendre(_Arg1, _Arg2, static_cast<double>(_Theta)));
	}

_NODISCARD double sph_neumann(unsigned, double);
_NODISCARD float sph_neumannf(unsigned, float);
_NODISCARD inline long double sph_neumannl(unsigned _Arg1, long double _Arg2)
	{
	return (_STD sph_neumann(_Arg1, static_cast<double>(_Arg2)));
	}

_NODISCARD double hypot(double _Dx, double _Dy, double _Dz) _NOEXCEPT;
_NODISCARD float hypot(float const _Dx, float const _Dy, float const _Dz) _NOEXCEPT;
_NODISCARD inline long double hypot(long double const _Dx, long double const _Dy,
	long double const _Dz) _NOEXCEPT
	{
	return (hypot(static_cast<double>(_Dx), static_cast<double>(_Dy), static_cast<double>(_Dz)));
	}

_STD_END
