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
#include <eve/constant/inf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/module/real/special/detail/bessel_kernel.hpp>
#include <eve/function/cyl_bessel_y0.hpp>
#include <eve/detail/hz_device.hpp>


namespace eve::detail
{

  template<real_scalar_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    if (is_eqz(nu)) return cyl_bessel_y0(x);
    std::cout << "scalarscalar" << std::endl;
    if(is_ltz(x)) return nan(as(x));
    if(is_eqz(x)) return minf(as(x));
   if  (sqr(x) < 10 * (inc(nu)))
    {
      std::cout << "scalarscalar small cyl_bessel_y" << std::endl;
      return kernel_ij_series(nu, x, T(1), 200);
    }
    else if (x > T(10000))
    {
      std::cout << "scalarscalar large cyl_bessel_y" << std::endl;
      if(x == inf(as(x))) return   zero(as(x));
      auto [j, y] = kernel_asymp_jy(nu, x);
      return y;
    }
    else
    {
      std::cout << "scalarscalar medium cyl_bessel_y" << std::endl;
      auto [j, jp, y, yp] = kernel_jy(nu, x);
      return y;
    }
  }

  template<real_scalar_value I, floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    std::cout << "scalarsimd" << std::endl;
     return cyl_bessel_y(T(nu), x);
  }

  template<real_simd_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
     std::cout << "simdscalar" << std::endl;
    using c_t = wide <T, cardinal_t<I>>;
    return cyl_bessel_y(convert(nu, as(x)), c_t(x));
  }

   template<integral_simd_value I, floating_real_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    std::cout << "integralany" << std::endl;
     return cyl_bessel_y(convert(nu, as(element_type_t<T>())), x);
  }

  template<floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), T nu, T x) noexcept
  {
    auto xnotinf   = x != inf(as<T>());
    if constexpr(has_native_abi_v<T>)
    {
      auto br_large = [xnotinf](auto nu, auto x){
        std::cout << "large cyl_bessel_y" << std::endl;
        x =  if_else(x > T(10000), x, T(10000));
        auto [j, y] = kernel_asymp_jy(nu, x);
        j = if_else(x == inf(as(x))
                   , zero
                   , j);
        return if_else(xnotinf, j, 0);
      };

      auto br_small = [](auto nu, auto x){
        std::cout << "small cyl_bessel_y" << std::endl;
        x =  if_else((sqr(x) < 10 * (inc(nu))), x, one);
        auto y = kernel_ij_series(nu, x, element_type_t<T>(1), 200);
        return if_else(is_eqz(x)
                      , minf(as(x))
                      , y);
      };

      auto br_medium = [](auto nu, auto x, auto notdone){
        std::cout << "medium cyl_bessel_y" << std::endl;
        x = if_else(notdone, x, one);
        auto [j, jp, y, yp] = kernel_jy(nu, x);
        return if_else(is_eqz(x)
                     , minf(as(x))
                     , y);
      };

      auto r = nan(as<T>()); //nan case treated here
      auto notdone = is_not_nan(x) && is_gez(x);
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
      return r;
    }
    else
      return apply_over(cyl_bessel_y, nu, x);
  }
}
