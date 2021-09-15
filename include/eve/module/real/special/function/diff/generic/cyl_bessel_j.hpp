//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/function/cyl_bessel_jy.hpp>
#include <eve/function/cyl_bessel_j.hpp>
#include <eve/function/derivative.hpp>


namespace eve::detail
{

  template<floating_real_value T, real_value N>
  EVE_FORCEINLINE constexpr auto cyl_bessel_j_(EVE_SUPPORTS(cpu_)
                                  , diff_type<1> const &
                                  , N const &nu
                                  , T const &x) noexcept
  {
    auto jm  = cyl_bessel_j(nu  , x);
    auto jm1 = cyl_bessel_j(nu-1, x);
    auto jpb =  jm-(nu/x)*jm1;
    return jpb;
  }
}
