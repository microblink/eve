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

#include <eve/as.hpp>
#include <eve/detail/implementation.hpp>
#include <eve/detail/is_native.hpp>
#include <eve/memory/aligned_ptr.hpp>

namespace eve::detail
{
  //================================================================================================
  // Common logical case
  //================================================================================================
  template<real_scalar_value T, typename N, native ABI>
  EVE_FORCEINLINE auto load(eve::as_<logical<wide<T, N, ABI>>> const &,
                            ABI const & mode,
                            logical<T> const* ptr) noexcept
  {
    using type = typename logical<wide<T, N, ABI>>::storage_type;
    auto raw_data = load(eve::as_<wide<T, N, ABI>> {}, mode, (T *)ptr);
    return type(raw_data);
  }

  template<typename T, typename N, std::size_t Align, native ABI>
  EVE_FORCEINLINE auto
  load(eve::as_<logical<wide<T, N, ABI>>> const &,
       ABI const &                          mode,
       aligned_ptr<logical<T> const, Align> ptr) noexcept
  {
    using type = typename logical<wide<T, N, ABI>>::storage_type;
    auto raw_data = load(eve::as_<wide<T, N, ABI>>{}, mode, aligned_ptr<T const, Align>((T const*)ptr.get()));
    return type(raw_data);
  }

  template<typename T, typename N, std::size_t Align, native ABI>
  EVE_FORCEINLINE auto
  load(eve::as_<logical<wide<T, N, ABI>>> const &,
       ABI const &                          mode,
       aligned_ptr<logical<T>, Align>       ptr) noexcept
  {
    using type = typename logical<wide<T, N, ABI>>::storage_type;
    auto raw_data = load(eve::as_<wide<T, N, ABI>>{}, mode, aligned_ptr<T, Align>((T *)ptr.get()));
    return type(raw_data);
  }
}
