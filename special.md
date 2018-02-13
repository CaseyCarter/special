## [Special math functions [sf.cmath]](http://eel.is/c++draft/sf.cmath)

* Figure out how to use bcp to extract boost.math

For most error cases, the functions throw `std::domain_error`. There are, however, cases that throw something else:
* `std::overflow_error`
* `std::underflow_error`
* `boost::math::rounding_error`
* `boost::math::evaluation_error`
I need to determine which if any of these need to be either remapped before propagating into user code, or if it's possible to use a custom Policy to do so.


### [Associated Laguerre polynomials [sf.cmath.assoc_laguerre]](http://eel.is/c++draft/sf.cmath.laguerre)

```c++
double assoc_laguerre(unsigned n, unsigned m, double x);
float assoc_laguerref(unsigned n, unsigned m, float x);
long double assoc_laguerrel(unsigned n, unsigned m, long double x);
```

[Boost: Laguerre (and Associated) Polynomials](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_poly/laguerre.html)

These are overloads (`T laguerre(unsigned, unsigned, T)`) in boost.

### [Associated Legendre functions [sf.cmath.assoc_legendre]](http://eel.is/c++draft/sf.cmath.legendre)

```c++
double assoc_legendre(unsigned l, unsigned m, double x);
float assoc_legendref(unsigned l, unsigned m, float x);
long double assoc_legendrel(unsigned l, unsigned m, long double x);
```

[Boost: Legendre (and Associated) Polynomials](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_poly/legendre.html)

These are overloads (`T legendre_p(unsigned, unsigned, T)`) in boost.

**WARNING: the boost functions include an additional (-1)^m factor that the standard functions do not.**

### [Beta function [sf.cmath.beta]](http://eel.is/c++draft/sf.cmath.beta)

```c++
double beta(double x, double y);
float betaf(float x, float y);
long double betal(long double x, long double y);
```

[Boost: Beta](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_beta/beta_function.html)

Again overloads `XX beta(T1, T2)`

### [Complete elliptic integral of the first kind [sf.cmath.comp_ellint_1]](http://eel.is/c++draft/sf.cmath.comp_ellint_1)

```c++
double comp_ellint_1(double k);
float comp_ellint_1f(float k);
long double comp_ellint_1l(long double k);
```

[Boost: Elliptic Integrals of the First Kind - Legendre Form](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/ellint/ellint_1.html)

Overloads `T ellint_1(T)`.

### [Complete elliptic integral of the second kind [sf.cmath.comp_ellint_2]](http://eel.is/c++draft/sf.cmath.comp_ellint_2)

```c++
double comp_ellint_2(double k);
float comp_ellint_2f(float k);
long double comp_ellint_2l(long double k);
```

[Boost: Elliptic Integrals of the Second Kind - Legendre Form](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/ellint/ellint_2.html)

Overloads `T ellint_2(T)`.

### [Complete elliptic integral of the third kind [sf.cmath.comp_ellint_3]](http://eel.is/c++draft/sf.cmath.comp_ellint_3)

```c++
double comp_ellint_3(double k, double nu);
float comp_ellint_3f(float k, float nu);
long double comp_ellint_3l(long double k, long double nu);
```

[Boost: Elliptic Integrals of the Third Kind - Legendre Form](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/ellint/ellint_3.html)

Overloads `XX ellint_3(T1, T2)`.

### [Regular modified cylindrical Bessel functions [sf.cmath.cyl_bessel_i]](http://eel.is/c++draft/sf.cmath.cyl_bessel_i)

```c++
double cyl_bessel_i(double nu, double x);
float cyl_bessel_if(float nu, float x);
long double cyl_bessel_il(long double nu, long double x);
```

[Boost: Modified Bessel Functions of the First and Second Kinds](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/bessel/mbessel.html)

Overloads: `XX cyl_bessel_i(T1, T2)`.

** WARNING: nervous about this definition.**

### [Cylindrical Bessel functions of the first kind [sf.cmath.cyl_bessel_j]](http://eel.is/c++draft/sf.cmath.cyl_bessel_j)

```c++
double cyl_bessel_j(double nu, double x);
float cyl_bessel_jf(float nu, float x);
long double cyl_bessel_jl(long double nu, long double x);
```

[Boost: Bessel Functions of the First and Second Kinds](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/bessel/bessel_first.html)

Overloads: `XX cyl_bessel_j(T1, T2)`.

### [Irregular modified cylindrical Bessel functions [sf.cmath.cyl_bessel_k]](http://eel.is/c++draft/sf.cmath.cyl_bessel_k)

```c++
double cyl_bessel_k(double nu, double x);
float cyl_bessel_kf(float nu, float x);
long double cyl_bessel_kl(long double nu, long double x);
```

[Boost: Modified Bessel Functions of the First and Second Kinds](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/bessel/mbessel.html)

Overloads: `XX cyl_bessel_k(T1, T2)`.

**NOTE:** The standard specifies special handling when the first parameter is integral, which would otherwise be undefined, which boost does not.

### [Cylindrical Neumann functions [sf.cmath.cyl_neumann]](http://eel.is/c++draft/sf.cmath.cyl_neumann)

```c++
double cyl_neumann(double nu, double x);
float cyl_neumannf(float nu, float x);
long double cyl_neumannl(long double nu, long double x);
```

[Boost: Bessel Functions of the First and Second Kinds](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/bessel/bessel_first.html)

Overloads: `XX cyl_neumann(T1, T2)`.

**NOTE:** The standard specifies special handling when the first parameter is integral, which would otherwise be undefined, which boost does not.

### [Incomplete elliptic integral of the first kind [sf.cmath.ellint_1]](http://eel.is/c++draft/sf.cmath.ellint_1)

```c++
double ellint_1(double k, double phi);
float ellint_1f(float k, float phi);
long double ellint_1l(long double k, long double phi);
```

[Boost: Elliptic Integrals of the First Kind - Legendre Form](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/ellint/ellint_1.html)

Overloads: `XX ellint_1(T1, T2)`.

### [Incomplete elliptic integral of the second kind [sf.cmath.ellint_2]](http://eel.is/c++draft/sf.cmath.ellint_2)

```c++
double ellint_2(double k, double phi);
float ellint_2f(float k, float phi);
long double ellint_2l(long double k, long double phi);
```

[Boost: Elliptic Integrals of the Second Kind - Legendre Form](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/ellint/ellint_2.html)

Overloads: `XX ellint_2(T1, T2)`.

### [Incomplete elliptic integral of the third kind [sf.cmath.ellint_3]](http://eel.is/c++draft/sf.cmath.ellint_3)

```c++
double ellint_3(double k, double nu, double phi);
float ellint_3f(float k, float nu, float phi);
long double ellint_3l(long double k, long double nu, long double phi);
```

[Boost: Elliptic Integrals of the Third Kind - Legendre Form](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/ellint/ellint_3.html)

Overloads: `XX ellint_3(T1, T2, T3)`.

### [Exponential integral [sf.cmath.expint]](http://eel.is/c++draft/sf.cmath.expint)

```c++
double expint(double x);
float expintf(float x);
long double expintl(long double x);
```

[Boost: Exponential Integral Ei](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/expint/expint_i.html)

Overloads: `T expint(T)`.

**NOTE:** the graph on the boost page indicates this is the correct function, despite that it is specified to calculate the *negative* of the correct function. I think it's a typo in the Boost documentation.

### [Hermite polynomials [sf.cmath.hermite]](http://eel.is/c++draft/sf.cmath.hermite)

```c++
double hermite(unsigned n, double x);
float hermitef(unsigned n, float x);
long double hermitel(unsigned n, long double x);
```

[Boost: Hermite Polynomials](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_poly/hermite.html)

Overloads: `T hermite(unsigned, T)`.

**NOTE:** Looks like there's another Boost documentation error: the formula says "d^2/dx^2" in contrast to literally all other sources that have "d^n/dx^n".

### [Laguerre polynomials [sf.cmath.laguerre]](http://eel.is/c++draft/sf.cmath.laguerre)

```c++
double laguerre(unsigned n, double x);
float laguerref(unsigned n, float x);
long double laguerrel(unsigned n, long double x);
```

[Boost: Laguerre (and Associated) Polynomials](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_poly/laguerre.html)

Overloads: `T laguerre(unsigned, T)`.

### [Legendre polynomials [sf.cmath.legendre]](http://eel.is/c++draft/sf.cmath.legendre)

```c++
double legendre(unsigned l, double x);
float legendref(unsigned l, float x);
long double legendrel(unsigned l, long double x);
```

[Boost: Legendre (and Associated) Polynomials](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_poly/legendre.html)

Overloads: `T legendre_p(int, T)`.

### [Riemann zeta function [sf.cmath.riemann_zeta]](http://eel.is/c++draft/sf.cmath.riemann_zeta)

```c++
double riemann_zeta(double x);
float riemann_zetaf(float x);
long double riemann_zetal(long double x);
```

[Boost: Riemann Zeta Function](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/zetas/zeta.html)

Overloads: `T zeta(T)`

**NOTE:** Verify behavior for x < 1.

### [Spherical Bessel functions of the first kind [sf.cmath.sph_bessel]](http://eel.is/c++draft/sf.cmath.sph_bessel)

```c++
double sph_legendre(unsigned l, unsigned m, double theta);
float sph_legendref(unsigned l, unsigned m, float theta);
long double sph_legendrel(unsigned l, unsigned m, long double theta);
```

[Boost: Spherical Harmonics](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_poly/sph_harm.html)

Overloads: `T spherical_harmonic_r(unsigned, int, T, 0)` (Note that `spherical_harmonic(l, m, theta, 0) == complex<decltype(theta)>(spherical_harmonic_r(l, m, theta, 0))`)

**NOTE:** The boost definitions DO NOT show the Condon-Shortley phase term (-1)^m because boost - unlike the Standard - includes that term in its definition of the associated legendre polynomials. The result is that THESE functions are equivalent to the standard C++ functions despite that the associated legendre polynomials are not.

**NOTE:** The boost functions accept int for the second parameter whereas the standard takes unsigned.

### [Spherical Neumann functions [sf.cmath.sph_neumann]](http://eel.is/c++draft/sf.cmath.sph_neumann)

```c++
double sph_neumann(unsigned n, double x);
float sph_neumannf(unsigned n, float x);
long double sph_neumannl(unsigned n, long double x);
```

[Boost: Spherical Bessel Functions of the First and Second Kinds](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/bessel/sph_bessel.html)

Overloads: `T sph_neumann(unsigned, T)` (I think there's another Boost documentation bug here: the boost function templates are depicted to take an additional non-deduced type parameter that is not used.)
