//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright 2019 Joel FALCOU
  Copyright 2019 Jean-Thierry LAPRESTE

  Licensed under the MIT License <http://opensource.org/licenses/MIT>.
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include <eve/function/is_lez.hpp>
#include <cstddef>

#define TYPE()        std::uint8_t
#define FUNCTION()    eve::is_lez
#define SAMPLES(N)    random<T>(N,0,200)

#include "bench.hpp"
