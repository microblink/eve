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
#include <eve/logical.hpp>
#include <eve/wide.hpp>

TTS_CASE_TPL( "Check enumerating constructor for arithmetic wide", EVE_TYPE )
{
  T ref([](auto i, auto) { return EVE_VALUE(i + 1); });

  T simd  = [&]<std::size_t... N>(std::index_sequence<N...>)
            {
              return T{(1+N)...};
            }( std::make_index_sequence<EVE_CARDINAL>());

  TTS_EQUAL(simd, ref);
}

TTS_CASE_TPL("Check enumerating constructor for wide of logical", EVE_TYPE)
{
  eve::logical<T> ref([](auto i, auto) { return eve::logical<EVE_VALUE>(i % 3 == 0); });
  eve::logical<T> simd  = [&]<std::size_t... N>(std::index_sequence<N...>)
                          {
                            return eve::logical<T>{(N%3 == 0)...};
                          }( std::make_index_sequence<EVE_CARDINAL>());

  TTS_EQUAL(simd, ref);
}