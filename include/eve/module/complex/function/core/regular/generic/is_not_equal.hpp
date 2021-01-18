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
#include <eve/function/bit_cast.hpp>
#include <eve/function/is_not_nan.hpp>
#include <eve/function/logical_or.hpp>
#include <eve/function/logical_and.hpp>
#include <eve/function/numeric.hpp>
#include <eve/logical.hpp>
#include <eve/traits/as_logical.hpp>
#include <type_traits>
#include <eve/concept/value.hpp>
#include <eve/concept/compatible.hpp>
#include <eve/detail/apply_over.hpp>
#include <eve/detail/skeleton_calls.hpp>

namespace eve::detail
{
  // -----------------------------------------------------------------------------------------------
  // regular case
  template<real_value T, real_value U>
  EVE_FORCEINLINE  auto is_not_equal_(EVE_SUPPORTS(cpu_)
                            , T const &a
                            , U const &b) noexcept
  requires compatible_values<T, U>
  {
    return arithmetic_call(is_not_equal, a, b);
  }

  template<real_scalar_value T>
  EVE_FORCEINLINE  auto is_not_equal_(EVE_SUPPORTS(cpu_)
                            , T const &a
                            , T const &b) noexcept
  {
    return as_logical_t<T>(a != b);
  }

  template<real_simd_value T>
  EVE_FORCEINLINE auto is_not_equal_(EVE_SUPPORTS(cpu_)
                            , T const &a
                            , T const &b) noexcept
  {
    return apply_over(is_not_equal, a, b);
  }

  template<real_value T, real_value U>
  EVE_FORCEINLINE  auto is_not_equal_(EVE_SUPPORTS(cpu_)
                            , logical<T> const &a
                            , logical<U> const &b) noexcept
  requires compatible_values<T, U>
  {
    return arithmetic_call(is_not_equal, a, b);
  }

  template<real_scalar_value T>
  EVE_FORCEINLINE  auto is_not_equal_(EVE_SUPPORTS(cpu_)
                            , logical<T> const &a
                            , logical<T> const &b) noexcept
  {
    return as_logical_t<T>(a.bits() != b.bits());
  }

  template<real_simd_value T>
  EVE_FORCEINLINE auto is_not_equal_(EVE_SUPPORTS(cpu_)
                            , logical<T> const &a
                            , logical<T> const &b) noexcept
  {
    return bit_cast(is_not_equal(a.bits(), b.bits()),as<logical<T>>());
  }
}

namespace eve
{
  template<value T, value U>
  EVE_FORCEINLINE auto operator !=(T const &v0, U const &v1) noexcept
  -> decltype( eve::is_not_equal(v0,v1) )
  requires compatible_values<T, U>
  {
    return eve::is_not_equal(v0, v1);
  }

  // This is needed to prevent clang using an internal operator of comparison over simd vector
  // types with different value types, as some architectures are not
  // aware of the accurate element type in the storage
  template<real_simd_value T, real_simd_value U>
  EVE_FORCEINLINE  auto operator != (T const &a
                                    , U const &b) noexcept
  requires different_value_type<T, U> = delete;
}