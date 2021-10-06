//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/detail/hz_device.hpp>
#include <eve/concept/value.hpp>
#include <eve/constant/pi.hpp>
#include <eve/constant/invpi.hpp>
#include <eve/function/all.hpp>
#include <eve/function/convert.hpp>
#include <eve/function/converter.hpp>
#include <eve/function/cospi.hpp>
#include <eve/function/cyl_bessel_y0.hpp>
#include <eve/function/cyl_bessel_y1.hpp>
#include <eve/function/factorial.hpp>
#include <eve/function/fma.hpp>
#include <eve/function/fms.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/log.hpp>
#include <eve/function/none.hpp>
#include <eve/function/pow.hpp>
#include <eve/function/rsqrt.hpp>
#include <eve/function/sincos.hpp>
#include <eve/function/sqr.hpp>
//#include <eve/module/real/special/detail/evaluate_rational.hpp>
#include <array>
#include <cassert>

namespace eve::detail
{

  template <real_scalar_value N, floating_real_scalar_value T>
  T bessel_yn_small_z(N n, T z)
  {
    //
    // See http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/06/01/04/01/02/
    //
    // Note that when called we assume that x < epsilon and n is a positive integer.
    //
    assert(n >= 0);
    assert((z < eps(as(z))));
    auto euler = Ieee_constant<T, 0x3f13c468U, 0x3fe2788cfc6fb619ULL>();

    if(n == 0)
    {
      return (2 * invpi(as(z))) * (eve::log(z*half(as(z))) + euler);
    }
    else if(n == 1)
    {
      return (z / pi(as(z))) * log(z / 2)
        - 2 / (pi(as(z)) * z)
        - (z / (2 * pi(as(z)))) * (1 - 2 * euler);
//      return (z* invpi(as(z))) * log(z*half(as(z))) - 2 * invpi(as(z)) * z - (z / (2 * pi(as(z))) * (1 - 2 * euler));                                                                              }
    }
    else if(n == 2)
    {
      return sqr(z/2) / (pi(as(z))) * log(z / 2)
        - (4 / (pi(as(z)) * z * z))
        - ((z * z) / (8 * pi(as(z)))) * (T(3)/2 - 2 * euler);
    }
    else
    {
      T p = pow(z / 2, n);
      T result = -((T(factorial(unsigned(n - 1))) / pi(as(z))));
      return result/p;
    }
  }

  template<scalar_value N, floating_real_scalar_value T>
  EVE_FORCEINLINE T cyl_bessel_yn_(EVE_SUPPORTS(cpu_), N nu, T x) noexcept
  {
    if constexpr(integral_scalar_value<N>) return cyl_bessel_yn(T(nu), x);
    else EVE_ASSERT(is_flint(nu), "nu is not a floating integral value");
    if (x < 0)  return nan(as(x));
    if (is_eqz(x)) return minf(as(x));
    if (x == inf(as(x))) return zero(as(x));
    if (nu < 0 ) return cospi(nu)*cyl_bessel_yn(-nu, x);
    if (x < eve::eps(as(x)))
    {
      return bessel_yn_small_z(nu, x);
    }
    else if (x > T(10000))//asymp
    {
      auto [j, y] = kernel_asymp_jy(T(nu), x);
      return y;
    }
    auto y0 = cyl_bessel_y0(x);
    if (nu == 0) return y0;
    auto y1 = cyl_bessel_y1(x);
    if (nu == 1) return y1;
    if (x == 0) return minf(as(x));
    // main case

    auto prev = y0;
    auto current = y1;
    int k = 1;
    T init = 2*rec(x);
    T mult = init;
    auto value = fms(mult, current, prev);
    prev = current;
    current = value;
    ++k;
    T factor = 1;
    if((mult > 1) && (eve::abs(current) > 1))
    {
      prev /= current;
      factor /= current;
      value /= current;
      current = 1;
    }
    while(k < nu)
    {
      mult = k*init;
      value = fms(mult, current, prev);
      prev = current;
      current = value;
      ++k;
    }
    return value/factor;
  }

  template<real_scalar_value I, floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_yn_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    return cyl_bessel_yn(T(nu), x);
  }

  template<real_simd_value I, floating_real_scalar_value T>
    EVE_FORCEINLINE auto cyl_bessel_yn_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    using c_t = wide <T, cardinal_t<I>>;
    return cyl_bessel_yn(convert(nu, as(x)), c_t(x));
  }

  template<real_simd_value I, floating_real_simd_value T>
  EVE_FORCEINLINE T cyl_bessel_yn_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
    requires(cardinal_v<T> == cardinal_v<I>)
  {
    using v_t =  eve::element_type_t<T>;
    if constexpr(integral_value<I>) return cyl_bessel_yn(convert(nu, as<v_t>()), x);
    else EVE_ASSERT(eve::all(is_flint(nu)), "some nu elements are not floating integral values");
    if constexpr(has_native_abi_v<T>)
    {

      auto br_small =  [](auto n, auto z)
        {
          auto euler = Ieee_constant<T, 0x3f13c468U, 0x3fe2788cfc6fb619ULL>();
          auto iseqzn = is_eqz(n);
          auto iseq1n = n == 1;
          auto iseq2n = n == 2;
          auto isgt2  = n > 2;
          auto r(zero(as(z)));
          if (eve::any(iseqzn))
          {
            r = (2 * invpi(as(z))) * (eve::log(z*half(as(z))) + euler);
          }
          if (eve::any(iseq1n))
          {
            auto r1 = (z / pi(as(z))) * log(z / 2) - 2 / (pi(as(z)) * z) - (z / (2 * pi(as(z)))) * (1 - 2 * euler);
            r = if_else(iseq1n, r1, r);
          }
          if(eve::any(iseq2n))
          {
            auto r2 = sqr(z/2) / (pi(as(z))) * log(z / 2) - (4 / (pi(as(z)) * z * z)) - ((z * z) / (8 * pi(as(z)))) * (T(3)/2 - 2 * euler);
            r = if_else(iseq1n, r2, r);
          }
          if(eve::any(isgt2))
          {
            T p = pow(z / 2, n);
            T result = -((convert(factorial(uint_(n - 1)), as<v_t>()) / pi(as(z))));
            r = if_else(isgt2, result/p, r);
          }
          return r;
        };

      auto br_large =  [](auto nu, auto x)
        {
          auto [j, y] = kernel_asymp_jy(nu, x);
          return y;
        };

      auto br_medium =  [](auto nu, auto x)
        {
          auto y0 = cyl_bessel_y0(x);
          if (eve::all(is_eqz(nu))) return y0;
          auto y1 = cyl_bessel_y1(x);
          if (eve::all(nu == 1)) return y1;
          // main case
          auto prev = y0;
          auto current = y1;
          int k = 1;
          T init = 2*rec(x);
          T mult = init;
          auto value = fms(mult, current, prev);
          prev = current;
          current = value;
          ++k;
          T factor(1);
          auto test = (mult > 1) && (eve::abs(current) > 1);
          if(eve::any(test))
          {
            current = if_else(test, rec(current)  , current);
            prev    = if_else(test, prev*current  , prev  );
            factor  = if_else(test, factor*current, factor);
            value   = if_else(test, value*current , value );
            current = if_else(test, one , current);
          }
          test = k < nu;
          while(eve::any(test))
          {
            mult  = if_else(test, k*init, mult);
            value = if_else(test, fms(mult, current, prev), value);
            prev = current;
            current = value;
            ++k;
            test = k < nu;
          }
          auto r = value/factor;
          r = if_else (is_eqz(nu), y0, r);
          r = if_else (nu == 1, y1, r);
          return r;

        };

      //reflection
      auto test = is_ltz(nu);
      if (eve::any(test))
      {
        auto f = if_else(test, cospi(nu), one);
        nu = if_else(test, -nu, nu);
        x  = if_else(test, -x, x);
        return f*cyl_bessel_yn(nu, x);
      }
      auto r = nan(as(x));
      test = is_eqz(x);
      r = if_else(test, minf(as(x)), r);
      auto notdone = is_gtz(x);
      if( eve::any(notdone) )
      {
        notdone = next_interval(br_small,  notdone, x <= eps(as(x)), r, nu, x);
        if( eve::any(notdone) )
        {
          notdone = next_interval(br_large,  notdone, x > T(10000), r, nu, x);
          if( eve::any(notdone) )
          {
            notdone = last_interval(br_medium,  notdone, r, nu, x);
          }
        }
      }
      r = if_else (x == inf(as(x)), zero, r);
      return r;
    }
    else return apply_over(cyl_bessel_yn, nu, x);
  }
}
