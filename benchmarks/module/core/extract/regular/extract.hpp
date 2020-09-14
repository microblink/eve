//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright 2020 Joel FALCOU
  Copyright 2020 Jean-Thierry LAPRESTE

  Licensed under the MIT License <http://opensource.org/licenses/MIT>.
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include <eve/function/extract.hpp>
#include <eve/constant/valmin.hpp>
#include <eve/constant/valmax.hpp>
#include <cmath>

int main()
{
  using I_VALUE  = eve::detail::as_integer_t<EVE_VALUE>;
  using I_TYPE   = eve::detail::as_integer_t<EVE_TYPE>;
  auto lmin = eve::valmin(eve::as<EVE_VALUE>());
  auto lmax = eve::valmax(eve::as<EVE_VALUE>());
  auto smin = I_VALUE(0);
  auto smax = I_VALUE(sizeof(EVE_VALUE)-1);

  auto arg0 = eve::bench::random_<EVE_VALUE>(lmin,lmax);
  auto arg1 = eve::bench::random_<I_VALUE>(smin,smax);

  eve::bench::experiment xp;
  run<eve::bench::types<EVE_TYPE,  I_VALUE>> (EVE_NAME(extract) , xp, eve::extract, arg0, arg1);


}
