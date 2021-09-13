//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/concept/value.hpp>
#include <eve/function/cyl_bessel_jn.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/is_ltz.hpp>
#include <eve/function/is_flint.hpp>
#include <eve/function/is_odd.hpp>

namespace eve::detail
{

  template<real_scalar_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    auto defined =  is_ltz(x) && eve::is_flint(nu);
    x = eve::abs(x);
    auto [j, jp, n, np] = cyl_bessel_jn(nu, x);
    j = defined && is_odd(nu) ? -x :x;
    return j;
  }

  template<real_scalar_value I, floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    return cyl_bessel_j(T(nu), x);
  }

  template<real_simd_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    using c_t = wide <T, cardinal_t<I>>;
    return cyl_bessel_j(convert(nu, as(x)), c_t(x));
  }

  template<real_value I, floating_real_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    auto defined =  is_ltz(x) && eve::is_flint(nu);
    x = if_else(defined, -x, x);
    auto [j, jp, n, np] = cyl_bessel_jn(nu, x);
    j = if_else(defined && is_odd(nu), -j, j);
    return j;
  }
}
