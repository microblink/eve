//==================================================================================================
/**
   EVE - Expressive Vector Engine
   Copyright 2020 Jean-Thierry LAPRESTE
   Copyright 2020 Joel FALCOU

   Licensed under the MIT License <http://opensource.org/licenses/MIT>.
   SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/concept/value.hpp>
#include <eve/constant/zero.hpp>
#include <eve/detail/apply_over.hpp>
#include <eve/detail/implementation.hpp>
#include <eve/function/trunc.hpp>

#include <tuple>

namespace eve::detail
{
  template<real_value T> EVE_FORCEINLINE constexpr auto modf_(EVE_SUPPORTS(cpu_), T a) noexcept
  {
    if constexpr( has_native_abi_v<T> )
    {
      if constexpr( floating_real_value<T> )
      {
        auto t = trunc(a);
        return std::make_tuple(a - t, t);
      }
      else
      {
        return std::make_tuple(zero(eve::as(a)), a);
      }
    }
    else
    {
      return apply_over2(modf, a);
    }
  }
}