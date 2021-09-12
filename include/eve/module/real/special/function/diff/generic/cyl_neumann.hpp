//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/function/cyl_bessel_jn.hpp>
#include <eve/function/derivative.hpp>

namespace eve::detail
{

  template<floating_real_value T, real_value N>
  EVE_FORCEINLINE constexpr T cyl_neumann_(EVE_SUPPORTS(cpu_)
                                          , diff_type<1> const &
                                          , N const &nu
                                          , T const &x) noexcept
  {
    auto [j, jp, n, np] = cyl_bessel_jn(nu, x);
    return np;
  }
}
