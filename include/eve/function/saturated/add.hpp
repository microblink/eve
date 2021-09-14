//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/function/add.hpp>
#include <eve/module/real/core/function/saturated/generic/add.hpp>

#if defined(EVE_INCLUDE_X86_HEADER)
#  include <eve/module/real/core/function/saturated/simd/x86/add.hpp>
#endif

namespace eve::tag { struct add_saturated_ {}; }

namespace eve::detail
{
  struct callable_add_saturated_;

  template<> struct decorate<tag::add_,saturated2_>
  {
    using type = callable_add_saturated_;
  };

  struct callable_add_saturated_  : make_callable<tag::add_saturated_   , callable_add_saturated_>
                                  , make_conditional<tag::add_saturated_, callable_add_saturated_>
  {
    static constexpr auto name = "saturated(add)";

    template<value... Tn>
    static EVE_FORCEINLINE constexpr auto call(Tn... an) noexcept
    {
      return add_(EVE_DISPATCH(), saturated_type{}, an...);
    }

    template<conditional_expr C, value... Tn>
    static EVE_FORCEINLINE constexpr auto conditional_call(C const& c, Tn... an) noexcept
    {
      return add_(EVE_DISPATCH(), c, saturated_type{}, an...);
    }
  };
}
