//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/concept/value.hpp>
#include <eve/function/all.hpp>
#include <eve/function/any.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/is_ltz.hpp>
#include <eve/function/is_flint.hpp>
#include <eve/function/is_not_nan.hpp>
#include <eve/function/is_odd.hpp>
#include <eve/function/nthroot.hpp>
#include <eve/constant/inf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/module/real/special/detail/bessel_kernel.hpp>
#include <eve/function/cyl_bessel_y0.hpp>
#include <eve/function/cyl_bessel_y1.hpp>
#include <eve/function/cyl_bessel_yn.hpp>
#include <eve/detail/hz_device.hpp>


namespace eve::detail
{

  template<real_scalar_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    std::cout << "bessel_jy scalar nu " << nu << " x " << x << std::endl;
    if (is_eqz(nu)) return cyl_bessel_y0(x);
    if (nu == I(1)) return cyl_bessel_y1(x);
    if (is_flint(nu)) return cyl_bessel_yn(nu, x);
    if(is_ltz(x)) return nan(as(x));
    if(is_eqz(x)) return minf(as(x));
//     if  (sqr(x) < 10 * (inc(nu)))
//     {
//       std::cout << "small " << std::endl;
//       return kernel_ij_series(T(nu), x, T(1), 200);
//     }
//    else
    if (x > T(10000))
    {
      std::cout << "big " << std::endl;
      if(x == inf(as(x))) return   zero(as(x));
      auto [j, y] = kernel_asymp_jy(nu, x);
      return y;
    }
    else
    {
     std::cout << "medium " << std::endl;
       auto [j, jp, y, yp] = kernel_jy(nu, x);
      return y;
    }
  }

  template<real_scalar_value I, floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
     return cyl_bessel_y(T(nu), x);
  }

  template<integral_simd_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    using c_t = wide <T, cardinal_t<I>>;
    return cyl_bessel_y(convert(nu, as(x)), c_t(x));
  }

  template<floating_real_simd_value U, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), U nu, T x) noexcept
  {
    return cyl_bessel_y(nu, U(x));
  }

  template<integral_simd_value I, floating_real_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    return cyl_bessel_yn(convert(nu, as(element_type_t<T>())), x);
  }

  template<floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_y_(EVE_SUPPORTS(cpu_), T nu, T x) noexcept
  {
    std::cout << "zut " << x << " " << nu << std::endl;
    if constexpr(has_native_abi_v<T>)
    {
      std::cout << "zut 1 " << std::endl;
      auto isflintnu = is_flint(nu);
      if (eve::all(isflintnu)) return cyl_bessel_yn(nu, x);
      std::cout << "zut 2 " << std::endl;
      auto xnotinf   = x != inf(as<T>());

      auto br_large = [xnotinf](auto nu, auto x){
        std::cout << "icitte large" << std::endl;
        auto xx =  if_else(x > T(1000), x, T(1000));
        auto [j, y] = kernel_asymp_jy(nu, xx);
 //        y = if_else(x == inf(as(x))
//                    , zero
//                    , y);
        return if_else(xnotinf, y, zero);
      };

//       auto br_small = [](auto nu, auto x){
//         x =  if_else((sqr(x) < 10 * (inc(nu))), x, one);
//         auto y = kernel_ij_series(nu, x, element_type_t<T>(1), 200);
//         return if_else(is_eqz(x)
//                       , minf(as(x))
//                       , y);
//       };

      auto br_medium = [](auto nu, auto x, auto notdone){
        std::cout << "medium " << x << std::endl;
        auto xx = if_else(notdone, x, one);
        auto [j, jp, y, yp] = kernel_jy(nu, xx);
        return if_else(is_eqz(x), minf(as(x)), y);
      };

      auto r = nan(as<T>()); //nan case treated here
      auto notdone = is_not_nan(x) && is_gez(x);
      std::cout << "1 " << x << " " << notdone << std::endl;
      if( eve::any(notdone) )
      {
//         notdone = next_interval(br_small,  notdone, (sqr(x) < 10 * (inc(nu))), r, nu, x);
//         if( eve::any(notdone) )
//         {
        notdone = next_interval(br_large,  notdone, x > T(1000), r, nu, x);
        std::cout << "2 " << x << " " << notdone << std::endl;
        if( eve::any(notdone) )
        {
          notdone = last_interval(br_medium,  notdone, r, nu, x, notdone);
        }
        //    }
      }
      auto tst = isflintnu && (x > T(1000));
      std::cout << "3 " << x << " " << tst << std::endl;
      if(eve::any(tst))
      {
        std::cout << "integral large" << std::endl;
        auto nnu = if_else(tst, nu, zero);
        auto xx  = if_else(tst, x,  nan(as(x)));
        return if_else(tst, cyl_bessel_yn(nnu, xx), r);
      }
      return r;
    }
    else
      return apply_over(cyl_bessel_y, nu, x);
  }
}
