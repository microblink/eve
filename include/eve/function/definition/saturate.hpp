//==================================================================================================
/**
  Copyright 2018 Joel FALCOU

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
**/
//==================================================================================================
#ifndef EVE_FUNCTION_DEFINITION_SATURATE_HPP_INCLUDED
#define EVE_FUNCTION_DEFINITION_SATURATE_HPP_INCLUDED

#include <eve/detail/overload.hpp>
#include <eve/detail/abi.hpp>
#include <eve/as.hpp>

namespace eve
{
  EVE_MAKE_CALLABLE(saturate_,saturate);

  template<typename Target, typename Arg>
  EVE_FORCEINLINE constexpr Arg saturate(Arg const& a0) noexcept
  {
    return saturate_(a0, as_<Target>{});
  }
}

#endif
