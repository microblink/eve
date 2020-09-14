//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright 2020 Joel FALCOU
  Copyright 2020 Jean-Thierry LAPRESTE

  Licensed under the MIT License <http://opensource.org/licenses/MIT>.
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include <cmath>
#include <eve/constant/valmin.hpp>
#include <eve/constant/valmax.hpp>
#include <eve/detail/meta/as_integer.hpp>

int main(int argc, char** argv)
{
  using I_TYPE =  eve::detail::as_integer_t<EVE_TYPE>;
  auto const std_ifrexp = [](auto x) { int e; auto m = std::ifrexp(x,  &e);
                                      return std::make_tuple(m, I_TYPE(e)); };
  auto l = eve::valmax(eve::as<EVE_TYPE>());
  EVE_REGISTER_BENCHMARK(std_ifrexp, EVE_TYPE
                        , eve::bench::random< EVE_TYPE>(-l, l));

  eve::bench::start_benchmarks(argc, argv);
}
