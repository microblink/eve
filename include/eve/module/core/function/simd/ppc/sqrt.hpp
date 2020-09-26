//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright 2020 Joel FALCOU
  Copyright 2020 Jean-Thierry LAPRESTE

  Licensed under the MIT License <http://opensource.org/licenses/MIT>.
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/detail/implementation.hpp>
#include <eve/function/rsqrt.hpp>
#include <eve/function/if_else.hpp>
#include <eve/constant/inf.hpp>
#include <type_traits>
#include <eve/concept/value.hpp>

namespace eve::detail
{
  template<floating_real_scalar_value T, typename N>
  EVE_FORCEINLINE wide<T, N, ppc_> sqrt_(EVE_SUPPORTS(vmx_), wide<T, N, ppc_> const &v0) noexcept
  {
    if constexpr(std::is_floating_point_v<T>)
    {
      if constexpr(current_api == vmx)
      {
        auto that = if_else(v0, v0 * rsqrt(v0), v0);
        if constexpr(platform::supports_invalids) { return if_else(v0 == inf(eve::as(v0)), v0, that); }
        else
        {
          return that;
        }
      }
      else
      {
        return vec_sqrt(v0.storage());
      }
    }
  }

  template<floating_real_scalar_value T, typename N>
  EVE_FORCEINLINE wide<T, N, ppc_>
                  sqrt_(EVE_SUPPORTS(vmx_), raw_type const &, wide<T, N, ppc_> const &v0) noexcept
  {
    if constexpr(std::is_floating_point_v<T>)
    {
      if constexpr(current_api == vmx) { return if_else(v0, v0 * raw(rsqrt)(v0), eve::zero); }
      else
      {
        return vec_sqrt(v0.storage());
      }
    }
  }
}

