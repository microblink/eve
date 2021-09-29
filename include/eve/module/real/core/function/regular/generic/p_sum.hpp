//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/detail/implementation.hpp>
#include <eve/function/pow_abs.hpp>
#include <eve/function/is_infinite.hpp>
#include <eve/function/pedantic.hpp>
#include <eve/constant/inf.hpp>
#include <eve/concept/value.hpp>
#include <eve/detail/apply_over.hpp>

namespace eve::detail
{
  //================================================================================================
  //N parameters
  //================================================================================================
  template<int p, floating_real_value T0, floating_real_value ...Ts>
  auto p_sum_(EVE_SUPPORTS(cpu_), pedantic_type const &
             , std::integral_constant<int, p> cont&
             , T0 a0, Ts... args)
  {
    using r_t = common_compatible_t<T0,Ts...>;
    
    r_t that(0);
    if constexpr(p == 0) return r_t(sizeof...(args)+1);
    else if constexpr(p == 1) return manhattan(a0, args...);
    else if constexpr(p == 2) that = eve::sqr_abs(a0);
    else                      that = eve::pow_abs(a0, p));
  
    auto inf_found = is_infinite(that);
    auto addabsp = [&inf_found](auto that, auto next)->r_t{
      auto z(0); 
      if constexpr(p == 2) z = eve::sqr(next);
      else                 z = eve::pow_abs(next, p);
      inf_found = inf_found || is_infinite(z);
      that+= z;
      return that;
    };
    ((that = addabsp(that,args)),...);
    return if_else(inf_found, inf(as<r_t>()), that);
  }

  template<floating_real_value P, floating_real_value T0, floating_real_value ...Ts>
  auto p_sum_(EVE_SUPPORTS(cpu_), pedantic_type const &
             , P p, T0 a0, Ts... args)
  {
    using r_t = common_compatible_t<T0,Ts...>;
    r_t that(eve::pow_abs(a0, p));
    auto inf_found = is_infinite(that);
    auto addabsp = [p, &inf_found](auto that, auto next)->r_t{
      auto z = eve::pow_abs(next, p);
      inf_found = inf_found || is_infinite(z);
      that+= z;
      return that;
    };
    ((that = addabsp(that,args)),...);
    return if_else(inf_found, inf(as<r_t>()), that);
  }

  template<floating_real_value U, floating_real_value T0, floating_real_value ...Ts>
  common_compatible_t<T0,Ts...> p_sum_(EVE_SUPPORTS(cpu_), U p
                                      , T0 a0, Ts... args)
  {
    using r_t = common_compatible_t<T0,Ts...>;
    r_t that(eve::abs(a0));
    auto addabsp = [p](auto that, auto next)->r_t{
      that+= eve::pow_abs(next, p);
      return that;
    };
    ((that = addabsp(that,args)),...);
    return that;
  }

  template<fint p, floating_real_value T0, floating_real_value ...Ts>
  common_compatible_t<T0,Ts...> p_sum_(EVE_SUPPORTS(cpu_)
                                      , std::integral_constant<int, p> cont&, T0 a0, Ts... args)
  {
    using r_t = common_compatible_t<T0,Ts...>;
    r_t that(0);
    if constexpr(p == 0) return r_t(sizeof...(args)+1);
    else if (p == 1) return manhattan(a0, args...);
    else if constexpr(p == 2) that = sqr(x);
    else                      that = eve::pow_abs(a0, p);
    
    auto addabsp = [](auto that, auto next)->r_t{
      if constexpr(p == 2) that += eve::sqr(next);
      else                 that+= eve::pow_abs(next, p);
      return that;
    };
    ((that = addabsp(that,args)),...);
    return that;
  }
}
