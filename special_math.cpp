#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/special_functions.hpp>
#include "special.hpp"

_STD_BEGIN
namespace {
template<class _Func> inline
	auto _Boost_call(const _Func _Fn)
	{
	_TRY_BEGIN
		return _Fn();
	_CATCH(boost::math::rounding_error&)
		_THROW(domain_error("FIXME: boost::math::rounding_error"));
	_CATCH(boost::math::evaluation_error&)
		_THROW(domain_error("FIXME: boost::math::evaluation_error"));
	_CATCH_END
	}
} // unnamed namespace

double assoc_laguerre(const unsigned _Pn, const unsigned _Pm, const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::laguerre(_Pn, _Pm, _Px); }));
	}

float assoc_laguerref(const unsigned _Pn, const unsigned _Pm, const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::laguerre(_Pn, _Pm, _Px); }));
	}

double assoc_legendre(const unsigned _Pl, const unsigned _Pm, const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::legendre_p(_Pl, _Pm, _Px); }));
	}

float assoc_legendref(const unsigned _Pl, const unsigned _Pm, const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::legendre_p(_Pl, _Pm, _Px); }));
	}

double beta(const double _Px, const double _Py)
	{
	return (_Boost_call([=]{ return boost::math::beta(_Px, _Py); }));
	}

float betaf(const float _Px, const float _Py)
	{
	return (_Boost_call([=]{ return boost::math::beta(_Px, _Py); }));
	}

double comp_ellint_1(const double _Pk)
	{
	return (_Boost_call([=]{ return boost::math::ellint_1(_Pk); }));
	}

float comp_ellint_1f(const float _Pk)
	{
	return (_Boost_call([=]{ return boost::math::ellint_1(_Pk); }));
	}

double comp_ellint_2(const double _Pk)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	return (_Boost_call([=]{ return boost::math::ellint_2(_Pk); }));
	}

float comp_ellint_2f(const float _Pk)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	return (_Boost_call([=]{ return boost::math::ellint_2(_Pk); }));
	}

double comp_ellint_3(const double _Pk, const double _Pnu)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	return (_Boost_call([=]{ return boost::math::ellint_3(_Pk, _Pnu); }));
	}

float comp_ellint_3f(const float _Pk, const float _Pnu)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	return (_Boost_call([=]{ return boost::math::ellint_3(_Pk, _Pnu); }));
	}

double cyl_bessel_i(const double _Pnu, const double _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_bessel_i(_Pnu, _Px); }));
}

float cyl_bessel_if(const float _Pnu, const float _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_bessel_i(_Pnu, _Px); }));
	}

double cyl_bessel_j(const double _Pnu, const double _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_bessel_j(_Pnu, _Px); }));
	}

float cyl_bessel_jf(const float _Pnu, const float _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_bessel_j(_Pnu, _Px); }));
	}

double cyl_bessel_k(const double _Pnu, const double _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_bessel_k(_Pnu, _Px); }));
	}

float cyl_bessel_kf(const float _Pnu, const float _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_bessel_k(_Pnu, _Px); }));
	}

double cyl_neumann(const double _Pnu, const double _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_neumann(_Pnu, _Px); }));
	}

float cyl_neumannf(const float _Pnu, const float _Px)
	{
	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::cyl_neumann(_Pnu, _Px); }));
	}

double ellint_1(const double _Pk, const double _Pphi)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pphi))
		{
		return (_Pphi);
		}

	return (_Boost_call([=]{ return boost::math::ellint_1(_Pk, _Pphi); }));
	}

float ellint_1f(const float _Pk, const float _Pphi)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pphi))
		{
		return (_Pphi);
		}

	return (_Boost_call([=]{ return boost::math::ellint_1(_Pk, _Pphi); }));
	}

double ellint_2(const double _Pk, const double _Pphi)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pphi))
		{
		return (_Pphi);
		}

	return (_Boost_call([=]{ return boost::math::ellint_2(_Pk, _Pphi); }));
	}

float ellint_2f(const float _Pk, const float _Pphi)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pphi))
		{
		return (_Pphi);
		}

	return (_Boost_call([=]{ return boost::math::ellint_2(_Pk, _Pphi); }));
	}

double ellint_3(const double _Pk, const double _Pnu, const double _Pphi)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Pphi))
		{
		return (_Pphi);
		}

	return (_Boost_call([=]{ return boost::math::ellint_3(_Pk, _Pnu, _Pphi); }));
	}

float ellint_3f(const float _Pk, const float _Pnu, const float _Pphi)
	{
	if (_CSTD isnan(_Pk))
		{
		return (_Pk);
		}

	if (_CSTD isnan(_Pnu))
		{
		return (_Pnu);
		}

	if (_CSTD isnan(_Pphi))
		{
		return (_Pphi);
		}

	return (_Boost_call([=]{ return boost::math::ellint_3(_Pk, _Pnu, _Pphi); }));
	}

double expint(const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::expint(_Px); }));
	}

float expintf(const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::expint(_Px); }));
	}

double hermite(const unsigned _Pn, const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::hermite(_Pn, _Px); }));
	}

float hermitef(const unsigned _Pn, const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::hermite(_Pn, _Px); }));
	}

double laguerre(const unsigned _Pn, const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::laguerre(_Pn, _Px); }));
	}

float laguerref(const unsigned _Pn, const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::laguerre(_Pn, _Px); }));
	}

double legendre(const unsigned _Pl, const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::legendre_p(_Pl, _Px); }));
	}

float legendref(const unsigned _Pl, const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::legendre_p(_Pl, _Px); }));
	}

double riemann_zeta(const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::zeta(_Px); }));
	}

float riemann_zetaf(const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::zeta(_Px); }));
	}

double sph_bessel(const unsigned _Pn, const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::sph_bessel(_Pn, _Px); }));
	}

float sph_besself(const unsigned _Pn, const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::sph_bessel(_Pn, _Px); }));
	}

double sph_legendre(const unsigned _Pl, const unsigned _Pm, const double _Ptheta)
	{
	if (_CSTD isnan(_Ptheta))
		{
		return (_Ptheta);
		}

	return (_Boost_call([=]{ return boost::math::spherical_harmonic_r(_Pl, _Pm, _Ptheta, 0.0); }));
	}

float sph_legendref(const unsigned _Pl, const unsigned _Pm, const float _Ptheta)
	{
	if (_CSTD isnan(_Ptheta))
		{
		return (_Ptheta);
		}

	return (_Boost_call([=]{ return boost::math::spherical_harmonic_r(_Pl, _Pm, _Ptheta, 0.0f); }));
	}

double sph_neumann(const unsigned _Pn, const double _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::sph_neumann(_Pn, _Px); }));
	}

float sph_neumannf(const unsigned _Pn, const float _Px)
	{
	if (_CSTD isnan(_Px))
		{
		return (_Px);
		}

	return (_Boost_call([=]{ return boost::math::sph_neumann(_Pn, _Px); }));
	}

namespace {
template<class _Ty> inline
	_Ty _Hypot3(_Ty _Dx, _Ty _Dy, _Ty _Dz)
	{
	static_assert(is_floating_point_v<_Ty>);
	_Dx = _CSTD fabs(_Dx);
	_Dy = _CSTD fabs(_Dy);
	_Dz = _CSTD fabs(_Dz);

	constexpr _Ty _Inf = numeric_limits<_Ty>::infinity();
	if (_Dx == _Inf || _Dy == _Inf || _Dz == _Inf)
		{
		return (_Inf);
		}

	if (_Dy > _Dx)
		{
		_STD swap(_Dx, _Dy);
		}

	if (_Dz > _Dx)
		{
		_STD swap(_Dx, _Dz);
		}

	constexpr _Ty _Eps = boost::math::tools::epsilon<_Ty>();
	if (_Dx * _Eps >= _Dy && _Dx * _Eps >= _Dz)
		{
		return (_Dx);
		}

	const auto _FracSq = [](const _Ty _Numerator, const _Ty _Denominator)
		{
		const _Ty result = _Numerator / _Denominator;
		return (result * result);
		};

	return (_Dx * _STD sqrt(1 + _FracSq(_Dy, _Dx) + _FracSq(_Dz, _Dx)));
	}
} // unnamed namespace

double hypot(const double _Dx, const double _Dy, const double _Dz)
	{
	return (_Hypot3<double>(_Dx, _Dy, _Dz));
	}

float hypot(const float _Dx, const float _Dy, const float _Dz)
	{
	return (_Hypot3<float>(_Dx, _Dy, _Dz));
	}
_STD_END
