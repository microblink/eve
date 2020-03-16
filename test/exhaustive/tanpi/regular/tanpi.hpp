//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright 2020 Joel FALCOU
  Copyright 2020 Jean-Thierry LAPRESTE

  Licensed under the MIT License <http://opensource.org/licenses/MIT>.
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include <eve/detail/function/tmp/boost_math_cospi.hpp>
#include <eve/detail/function/tmp/boost_math_sinpi.hpp>
#include <eve/function/tanpi.hpp>
#include <eve/constant/invpi.hpp>
#include <eve/constant/valmin.hpp>
#include <eve/constant/valmax.hpp>    
#include <eve/constant/maxflint.hpp>
#include <eve/constant/nan.hpp>
#include <eve/function/is_flint.hpp>
#include <eve/function/frac.hpp>
#include <tts/tests/range.hpp>
#include "measures.hpp"
#include "producers.hpp"
#include <cmath>

TTS_CASE("wide exhaustive check on tanpi")
{
  auto my_stdtanpi =  tts::vectorize<EVE_TYPE>([](auto x){return (eve::frac(std::abs(x)) == EVE_VALUE(0.5))
                                               ?  eve::Nan<EVE_VALUE>()
                                               : boost::math::sin_pi(x)/boost::math::cos_pi(x); });
  eve::exhaustive_producer<EVE_TYPE> p(eve::Valmin<EVE_VALUE>(), eve::Valmax<EVE_VALUE>());
  TTS_RANGE_CHECK(p, my_stdtanpi, eve::tanpi); 
}
