//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/detail/hz_device.hpp>
#include <eve/concept/value.hpp>
#include <eve/function/fma.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/none.hpp>
#include <eve/function/rec.hpp>

namespace eve::detail
{

  template <class T, class U, class V>
  EVE_FORCEINLINE V evaluate_rational(const T num, const U den, V z) noexcept
  {
    auto N =  num.size();
    auto eval_small = [&num, &den, N](auto z)
      {
        V s1(num[N-1]);
        V s2(den[N-1]);
        for(int i = (int)N - 2; i >= 0; --i)
        {
          s1 = fma(s1, z, num[i]);
          s2 = fma(s2, z, den[i]);
        }
        return s1/s2;
      };
    auto eval_large = [&num, &den, N](auto z)
      {
        z = rec(z);
        V s1(num[0]);
        V s2(den[0]);
        for(unsigned i = 1; i < N; ++i)
        {
          s1 = fma(s1, z, num[i]);
          s2 = fma(s2, z, den[i]);
        }
        return s1/s2;
       };
    auto test = z <= V(1);
    if(eve::all(test)) return eval_small(z);
    else if(eve::none(test)) return  eval_large(z);
    else return if_else(test,  eval_small(z),  eval_large(z));
  }
}
