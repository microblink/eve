//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/concept/value.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/is_ltz.hpp>
#include <eve/function/is_flint.hpp>
#include <eve/function/is_not_nan.hpp>
#include <eve/function/is_odd.hpp>
#include <eve/module/real/special/detail/bessel_kernel.hpp>
#include <eve/detail/hz_device.hpp>
#include <iostream>

namespace eve::detail
{

  template<real_scalar_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    std::cout << "scalarscalar" << std::endl;
    if(x == eve::inf(as(x))) return zero(as(x));
    auto defined =  is_ltz(x) && eve::is_flint(nu);
    x = eve::abs(x);
    if  (sqr(x) < 10 * (inc(nu)))
    {
      std::cout << "scalarscalar small" << std::endl;
      return kernel_ij_series(nu, x, T(-1), 200);
    }
    else if (x > T(10000))
    {
      std::cout << "scalarscalar large" << std::endl;
      if(x == eve::inf(as(x))) return  defined ? zero(as(x)): nan(as(x));
      auto [j, y] = kernel_asymp_jy(nu, x);
      return j;
    }
    else
    {
      std::cout << "scalarscalar medium" << std::endl;
      auto [j, jp, n, np] = kernel_jy(nu, x);
      j = defined && is_odd(nu) ? -j :j;
      return j;
    }
  }

  template<real_scalar_value I, floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    std::cout << "scalarsimd" << std::endl;
     return cyl_bessel_j(T(nu), x);
  }

  template<real_simd_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
     std::cout << "simdscalar" << std::endl;
    using c_t = wide <T, cardinal_t<I>>;
    return cyl_bessel_j(convert(nu, as(x)), c_t(x));
  }

   template<integral_simd_value I, floating_real_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    std::cout << "integralany" << std::endl;
     return cyl_bessel_j(convert(nu, as(element_type_t<T>())), x);
  }

  template<floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), T nu, T x) noexcept
  {
    auto xneg    = is_ltz(x);
    auto defined = xneg && eve::is_flint(nu);
    if constexpr(has_native_abi_v<T>)
    {
      auto br_large = [xneg, defined](auto nu, auto x){
        std::cout << "large" << std::endl;
        x =  if_else(x > T(10000), x, T(10000));
        auto [j, y] = kernel_asymp_jy(nu, x);
        j = if_else(x == eve::inf(as(x))
                   , if_else( xneg
                            , if_else(defined
                                     , zero
                                     , nan(as(x)))
                            , zero)
                   , j);
        return j;
      };

      auto br_small = [](auto nu, auto x){
       std::cout << "small" << std::endl;
         x =  if_else((sqr(x) < 10 * (inc(nu))), x, one);
       auto j = kernel_ij_series(nu, x, element_type_t<T>(-1), 200);
       return if_else(is_eqz(x)
                     , if_else(is_eqz(nu)
                              , one(as(x))
                              , zero)
                     , j);
      };

      auto br_medium = [](auto nu, auto x, auto notdone){
        std::cout << "medium" << std::endl;
       x = if_else(notdone, x, zero);
        auto [j, jp, n, np] = kernel_jy(nu, x);
        return j;
      };

      x = if_else(defined, -x, x);
      auto r = nan(as<T>()); //nan case treated here
      auto notdone = is_not_nan(x);
      if( eve::any(notdone) )
      {
        notdone = next_interval(br_small,  notdone, (sqr(x) < 10 * (inc(nu))), r, nu, x);
        if( eve::any(notdone) )
        {
          notdone = next_interval(br_large,  notdone, x > T(10000), r, nu, x);
          if( eve::any(notdone) )
          {
            notdone = last_interval(br_medium,  notdone, r, nu, x, notdone);
          }
        }
      }
      r = if_else(defined && is_odd(nu), -r, r);

      return r;
    }
    else
      return apply_over(cyl_bessel_j, nu, x);
  }
}
