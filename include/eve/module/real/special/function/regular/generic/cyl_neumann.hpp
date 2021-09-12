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


namespace eve::detail
{

  template<real_value I, floating_real_value T>
  EVE_FORCEINLINE auto cyl_neumann_(EVE_SUPPORTS(cpu_), I nu, T x) noexcept
  {
    auto [j, jp, n, np] = cyl_bessel_jn(nu, x);
    return n;
  }
}
