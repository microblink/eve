//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright 2018 Joel FALCOU

  Licensed under the MIT License <http://opensource.org/licenses/MIT>.
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#ifndef EVE_MEMORY_ALIGN_HPP_INCLUDED
#define EVE_MEMORY_ALIGN_HPP_INCLUDED

#include <eve/memory/power_of_2.hpp>
#include <type_traits>
#include <cstdint>

namespace eve
{
  struct over   { std::size_t value; };
  struct under  { std::size_t value; };

  template< typename T
          , typename = std::enable_if_t<std::is_integral_v<T>>
          >
  constexpr auto align( T value, over alignment) noexcept
  {
    assert(is_power_of_2(alignment.value));
    return (value+alignment.value-1) & ~(alignment.value-1);
  }

  template< typename T
          , typename = std::enable_if_t<std::is_integral_v<T>>
          >
  constexpr auto align( T value, under alignment) noexcept
  {
    assert(is_power_of_2(alignment.value));
    return value & ~(alignment.value-1);
  }

  template<typename T> constexpr auto align( T* ptr, over alignment) noexcept
  {
    return reinterpret_cast<T*>( align(reinterpret_cast<std::uintptr_t>(ptr),alignment) );
  }

  template<typename T> constexpr auto align( T* ptr, under alignment) noexcept
  {
    return reinterpret_cast<T*>( align(reinterpret_cast<std::uintptr_t>(ptr),alignment) );
  }
}

#endif