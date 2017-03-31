
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#ifndef _XSPECIAL_MATH_H
#define _XSPECIAL_MATH_H
#ifndef RC_INVOKED
#include <yvals.h>

//#include <boost/math/special_functions/math_fwd.hpp>
//#include <boost/math/special_functions/factorials.hpp>
//#include <boost/math/tools/config.hpp>

 #pragma pack(push,_CRT_PACKING)
 #pragma warning(push,_STL_WARNING_LEVEL)
 #pragma warning(disable: _STL_DISABLED_WARNINGS)
 #pragma push_macro("new")
 #undef new

_STD_BEGIN

template<class _Ty>
constexpr _Ty _Legendre_next(unsigned l, _Ty x, _Ty Pl, _Ty Plm1)
	{	// Recurrence relation for legendre P and Q polynomials
	return (((2 * l + 1) * x * Pl - l * Plm1) / (l + 1));
	}

template<class _Ty> inline
_Ty _Legendre_imp(unsigned l, _Ty x)
	{	// implement Legendre P polynomials via recurrence
	if (x < -1 || x > 1)
		{
		_THROW(domain_error, "The Legendre Polynomial is only defined for -1 <= x <= 1");
		}

	// A solution of the first kind (P):
	_Ty p0 = 1, p1 = x;
	if (l == 0)
		{
		return (p0);
		}

	for (unsigned n = 1; n < l; ++n)
		{
		_STD swap(p0, p1);
		p1 = _Legendre_next(n, x, p0, p1);
		}

	return (p1);
	}

inline double legendre(unsigned l, double x)
	{
	return _Legendre_imp(l, x);
	}

inline float legendref(unsigned l, float x)
	{
	return _Legendre_imp(l, x);
	}

inline long double legendrel(unsigned l, long double x)
	{
	return _Legendre_imp(l, x);
	}

template<class _Ty>
constexpr _Ty _Legendre_next(unsigned l, unsigned m, _Ty x, _Ty Pl, _Ty Plm1)
	{	// recurrence for associated polynomials
	return (((2 * l + 1) * x * Pl - (l + m) * Plm1) / (l + 1 - m));
	}

template<class _Ty> inline
_Ty _Legendre_p_imp(int l, int m, _Ty x, _Ty sin_theta_power)
	{	// Legendre P associated polynomial
	// Error handling:
	if ((x < -1) || (x > 1))
		return policies::raise_domain_error<_Ty>(
		"boost::math::legendre_p<%1%>(int, int, %1%)",
			"The associated Legendre Polynomial is defined for"
			" -1 <= x <= 1, but got x = %1%.", x, pol);
	// Handle negative arguments first:
	if (l < 0)
		return _Legendre_p_imp(-l-1, m, x, sin_theta_power, pol);
	if (m < 0)
	{
		int sign = (m&1) ? -1 : 1;
		return sign * boost::math::tgamma_ratio(static_cast<_Ty>(l+m+1), static_cast<_Ty>(l+1-m), pol) * _Legendre_p_imp(l, -m, x, sin_theta_power, pol);
	}
	// Special cases:
	if (m > l)
		return 0;
	if (m == 0)
		return boost::math::legendre_p(l, x, pol);

	_Ty p0 = boost::math::double_factorial<_Ty>(2 * m - 1, pol) * sin_theta_power;

	if (m&1)
		p0 *= -1;
	if (m == l)
		return p0;

	_Ty p1 = x * (2 * m + 1) * p0;

	int n = m + 1;

	while(n < l)
	{
		_STD swap(p0, p1);
		p1 = _Legendre_next(n, m, x, p0, p1);
		++n;
	}
	return p1;
}

template<class _Ty, class Policy> inline
_Ty _Legendre_p_imp(int l, int m, _Ty x, const Policy& pol)
{
	BOOST_MATH_STD_USING
	// TODO: we really could use that mythical "pow1p" function here:
	return _Legendre_p_imp(l, m, x, static_cast<_Ty>(pow(1 - x*x, _Ty(abs(m))/2)), pol);
}

}

template<class _Ty, class Policy> inline
typename tools::promote_args<_Ty>::type legendre_p(int l, int m, _Ty x, const Policy& pol)
{
	typedef typename tools::promote_args<_Ty>::type result_type;
	typedef typename policies::evaluation<result_type, Policy>::type value_type;
	return policies::checked_narrowing_cast<result_type, Policy>(detail::_Legendre_p_imp(l, m, static_cast<value_type>(x), pol), "bost::math::legendre_p<%1%>(int, int, %1%)");
}

template<class _Ty> inline
typename tools::promote_args<_Ty>::type legendre_p(int l, int m, _Ty x)
{
	return boost::math::legendre_p(l, m, x, policies::policy<>());
}

_STD_END

#endif // _XSPECIAL_MATH_H
