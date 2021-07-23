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
#include <eve/detail/apply_over.hpp>
#include <eve/detail/skeleton_calls.hpp>
#include <eve/concept/value.hpp>
#include <eve/concept/compatible.hpp>
#include <eve/function/bit_cast.hpp>
#include <eve/traits/as_logical.hpp>

namespace eve::detail
{
  template<int k, integral_value T, value U, value V>
  EVE_FORCEINLINE auto bit_ternary_(EVE_SUPPORTS(cpu_)
                                   , T const &a
                                   , U const &b
                                   , V const &c
                                   , std::integral_constant<int, k> const & ik) noexcept
  requires bit_compatible_values<T, U> && bit_compatible_values<T, V> && bit_compatible_values<U, V>
  {
    using vt_t = element_type_t<T>;
    if constexpr((sizeof(T) == sizeof(U)) && (sizeof(U) == sizeof(V)))    //all three are of the same bit size
    {
      if constexpr(has_native_abi_v<T> && has_native_abi_v<U> && has_native_abi_v<V>)
        return bit_ternary(a, bit_cast(b, as(a)), bit_cast(c, as(b)), ik);
      else
      {
        auto f =  [ik](auto a,  auto b,  auto c){ return bit_ternary(a, b, c, ik); };
        if constexpr(has_aggregated_abi_v<T>||has_aggregated_abi_v<U>||has_aggregated_abi_v<V>)  return aggregate(f, a, b, c);
        else if constexpr(has_emulated_abi_v<T>||has_emulated_abi_v<U>||has_emulated_abi_v<V>)   return map(f, a, b, c);
      }
    }
    else if constexpr(scalar_value<T> && scalar_value<U> && simd_value<V>)  //T, U are scalar so of the same bit size, V is simd
    {
      using r_t = as_wide_t<T, cardinal_t<V>>;
      auto  aa = r_t(bit_cast(a, as<vt_t>()));
      auto  bb = r_t(bit_cast(b, as<vt_t>()));
      auto  cc =   bit_cast(c, as<r_t>());
      return bit_ternary(aa, bb, cc, ik);
    }
    else if constexpr(scalar_value<T> && scalar_value<V> && simd_value<U>)  //T, V are scalar so of the same bit size, U is simd
    {
      using r_t = as_wide_t<T, cardinal_t<U>>;
      auto  aa = r_t(bit_cast(a, as<vt_t>()));
      auto  cc = r_t(bit_cast(c, as<vt_t>()));
      return bit_ternary(aa, b, cc, ik);
    }
    else if constexpr(scalar_value<U> && scalar_value<V> && simd_value<T>) //U, V are scalar so of the same bit size, T is simd
    {
      using r_t = T;
      auto  aa = bit_cast(a, as<r_t>());
      auto  bb = r_t(bit_cast(b, as<vt_t>()));
      auto  cc = r_t(bit_cast(b, as<vt_t>()));
      return bit_ternary(aa, bb, cc, ik);
    }
    else if constexpr(simd_value<U> && simd_value<V> && scalar_value<T>) //U, V are simd so of the same bit size, T is scalar
    {
      using r_t = as_wide_t<T, cardinal_t<U>>;
      auto  aa = r_t(bit_cast(a, as<vt_t>()));
      auto  cc = bit_cast(c, as<r_t>());
      return bit_ternary(aa, b, c, ik);
    }
    else if constexpr(simd_value<T> && simd_value<U> && scalar_value<V>) //U, T are simd so of the same bit size, V is scalar
    {
      using r_t = T;
      auto  bb = bit_cast(b, as<r_t>());
      auto  cc = r_t(bit_cast(c, as<vt_t>()));
      return bit_ternary(a, bb, cc, ik);
    }
    else if constexpr(simd_value<T> && simd_value<V> && scalar_value<U>) //T, V are simd so of the same bit size, U is scalar
    {
      using r_t = T;
      auto  bb = r_t(b);
      auto  cc = bit_cast(c, as<r_t>());
      return bit_ternary(a, bb, cc, ik);
    }
//     else if constexpr(simd_value<T> && simd_value<U> && simd_value<V>) // all are simd so of the same bit size
//     {
//       if constexpr(has_native_abi_v<T> && has_native_abi_v<U> && has_native_abi_v<V>)
//       {
//         using r_t = U;
//         auto  aa = bit_cast(a, as<r_t>());
//         auto  cc = bit_cast(c, as<r_t>());
//         return bit_ternary(aa, b, cc, ik); // generally already taken by arch specific intrisics
//       }
//       else return apply_over(op, a, b, c);
//     }
  }


  // this is adapted from Samuel neves ternary logic for sse avx etc.
  template < int k, integral_value T>
  EVE_FORCEINLINE T bit_ternary_(EVE_SUPPORTS(cpu_)
                                , [[maybe_unused]] T const & a
                                , [[maybe_unused]] T const & b
                                , [[maybe_unused]] T const & c
                                , std::integral_constant<int, k> const &) noexcept
  requires(k < 256)
  {
    if constexpr(k == 0x00)  // function=0, lowered=0, set=intel
    {
      return Zero(as(a));
    }
    else if constexpr(k == 0x01)  // function=(a nor (b or c)), lowered=((a or (b or c)) xor 1), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_or(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k == 0x02) // function=(c and (b nor a)), lowered=((b or a) notand c), set=intel
    {
      const T t0 = bit_or(b, a);
      const T t1 = bit_notand(t0, c);
      return t1;
    }
    else if constexpr(k == 0x03) // function=(b nor a), lowered=((b or a) xor 1), set=intel
    {
      const T t0 = bit_or(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k == 0x04) // function=(b and (a nor c)), lowered=((a or c) notand b), set=intel
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_notand(t0, b);
      return t1;
    }
    else if constexpr(k == 0x05) // function=(c nor a), lowered=((c or a) xor 1), set=intel
    {
      const T t0 = bit_or(c, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0x06) // function=(not (a) and (b xor c)), lowered=(a notand (b xor c)), set=automat
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_notand(a, t0);
      return t1;
    }
    else if constexpr(k ==  0x07) // function=(a nor (b and c)), lowered=((a or (b and c)) xor 1), set=intel
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_or(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x08) // function=((not (a) and b) and c), lowered=((a notand b) and c), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T t1 = bit_and(t0, c);
      return t1;
    }
    else if constexpr(k ==  0x09) // function=(a nor (b xor c)), lowered=((a or (b xor c)) xor 1), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_or(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x0a) // function=(c and not (a)), lowered=(a notand c), set=intel
    {
      const T t0 = bit_notand(a, c);
      return t0;
    }
    else if constexpr(k ==  0x0b) // function=(not (a) and ((b xor 1) or c)), lowered=(a notand ((b xor 1) or c)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      const T t1 = bit_or(t0, c);
      const T t2 = bit_notand(a, t1);
      return t2;
    }
    else if constexpr(k ==  0x0c) // function=(b and not (a)), lowered=(a notand b), set=intel
    {
      const T t0 = bit_notand(a, b);
      return t0;
    }
    else if constexpr(k ==  0x0d) // function=(not (a) and (b or (c xor 1))), lowered=(a notand (b or (c xor 1))), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      const T t1 = bit_or(b, t0);
      const T t2 = bit_notand(a, t1);
      return t2;
    }
    else if constexpr(k ==  0x0e) // function=(not (a) and (b or c)), lowered=(a notand (b or c)), set=automat
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_notand(a, t0);
      return t1;
    }
    else if constexpr(k ==  0x0f) // function=not (a), lowered=(a xor 1), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      return t0;
    }
    else if constexpr(k ==  0x10) // function=(a and (b nor c)), lowered=((b or c) notand a), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_notand(t0, a);
      return t1;
    }
    else if constexpr(k ==  0x11) // function=(c nor b), lowered=((c or b) xor 1), set=intel
    {
      const T t0 = bit_or(c, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0x12) // function=(not (b) and (a xor c)), lowered=(b notand (a xor c)), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_notand(b, t0);
      return t1;
    }
    else if constexpr(k ==  0x13) // function=(b nor (a and c)), lowered=((b or (a and c)) xor 1), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_or(b, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x14) // function=(not (c) and (a xor b)), lowered=(c notand (a xor b)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T t1 = bit_notand(c, t0);
      return t1;
    }
    else if constexpr(k ==  0x15) // function=(c nor (b and a)), lowered=((c or (b and a)) xor 1), set=intel
    {
      const T t0 = bit_and(b, a);
      const T t1 = bit_or(c, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x16) // function=(a ? (b nor c) : (b xor c)), lowered=(((b or c) notand a) or (a notand (b xor c))), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_notand(t0, a);
      const T t2 = bit_xor(b, c);
      const T t3 = bit_notand(a, t2);
      const T t4 = bit_or(t1, t3);
      return t4;
    }
    else if constexpr(k ==  0x17) // function=((b nor c) or (not (a) and (b xor c))), lowered=(((b or c) xor 1) or (a notand (b xor c))), set=optimized
    {
      const T t0 = bit_or(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_xor(b, c);
      const T t3 = bit_notand(a, t2);
      const T t4 = bit_or(t1, t3);
      return t4;
    }
    else if constexpr(k ==  0x18) // function=((a xor b) and (a xor c)), lowered=((a xor b) and (a xor c)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x19) // function=(not ((a and b)) and (b xor (c xor 1))), lowered=((a and b) notand (b xor (c xor 1))), set=automat
    {
      const T t0 = bit_and(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(b, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x1a) // function=(not ((a and b)) and (a xor c)), lowered=((a and b) notand (a xor c)), set=automat
    {
      const T t0 = bit_and(a, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x1b) // function=(c ? not (a) : not (b)), lowered=((a notand c) or (c notand (b xor 1))), set=intel
    {
      const T t0 = bit_notand(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_notand(c, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x1c) // function=(not ((a and c)) and (a xor b)), lowered=((a and c) notand (a xor b)), set=automat
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x1d) // function=(b ? not (a) : not (c)), lowered=((a notand b) or (b notand (c xor 1))), set=intel
    {
      const T t0 = bit_notand(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_notand(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x1e) // function=(a xor (b or c)), lowered=(a xor (b or c)), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_xor(a, t0);
      return t1;
    }
    else if constexpr(k ==  0x1f) // function=(a nand (b or c)), lowered=((a and (b or c)) xor 1), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_and(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x20) // function=((not (b) and a) and c), lowered=((b notand a) and c), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T t1 = bit_and(t0, c);
      return t1;
    }
    else if constexpr(k ==  0x21) // function=(b nor (a xor c)), lowered=((b or (a xor c)) xor 1), set=intel
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_or(b, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x22) // function=(c and not (b)), lowered=(b notand c), set=intel
    {
      const T t0 = bit_notand(b, c);
      return t0;
    }
    else if constexpr(k ==  0x23) // function=(not (b) and ((a xor 1) or c)), lowered=(b notand ((a xor 1) or c)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      const T t1 = bit_or(t0, c);
      const T t2 = bit_notand(b, t1);
      return t2;
    }
    else if constexpr(k ==  0x24) // function=((a xor b) and (b xor c)), lowered=((a xor b) and (b xor c)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x25) // function=(not ((a and b)) and (a xor (c xor 1))), lowered=((a and b) notand (a xor (c xor 1))), set=automat
    {
      const T t0 = bit_and(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x26) // function=(not ((a and b)) and (b xor c)), lowered=((a and b) notand (b xor c)), set=automat
    {
      const T t0 = bit_and(a, b);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x27) // function=(c ? not (b) : not (a)), lowered=((b notand c) or (c notand (a xor 1))), set=intel
    {
      const T t0 = bit_notand(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_notand(c, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x28) // function=(c and (b xor a)), lowered=(c and (b xor a)), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T t1 = bit_and(c, t0);
      return t1;
    }
    else if constexpr(k ==  0x29) // function=(c ? (b xor a) : (b nor a)), lowered=((c and (b xor a)) or (c notand ((b or a) xor 1))), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T t1 = bit_and(c, t0);
      const T t2 = bit_or(b, a);
      const T c1 = allbits(as<T>());
      const T t3 = bit_xor(t2, c1);
      const T t4 = bit_notand(c, t3);
      const T t5 = bit_or(t1, t4);
      return t5;
    }
    else if constexpr(k ==  0x2a) // function=(c and (b nand a)), lowered=((b and a) notand c), set=intel
    {
      const T t0 = bit_and(b, a);
      const T t1 = bit_notand(t0, c);
      return t1;
    }
    else if constexpr(k ==  0x2b) // function=(c ? (b nand a) : (b nor a)), lowered=(((b and a) notand c) or (c notand ((b or a) xor 1))), set=intel
    {
      const T t0 = bit_and(b, a);
      const T t1 = bit_notand(t0, c);
      const T t2 = bit_or(b, a);
      const T c1 = allbits(as<T>());
      const T t3 = bit_xor(t2, c1);
      const T t4 = bit_notand(c, t3);
      const T t5 = bit_or(t1, t4);
      return t5;
    }
    else if constexpr(k ==  0x2c) // function=((b or c) and (a xor b)), lowered=((b or c) and (a xor b)), set=automat
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x2d) // function=(a xor (b or (c xor 1))), lowered=(a xor (b or (c xor 1))), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      const T t1 = bit_or(b, t0);
      const T t2 = bit_xor(a, t1);
      return t2;
    }
    else if constexpr(k ==  0x2e) // function=((b or c) xor (a and b)), lowered=((b or c) xor (a and b)), set=optimized
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_and(a, b);
      const T t2 = bit_xor(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x2f) // function=((not (b) and c) or (a xor 1)), lowered=((b notand c) or (a xor 1)), set=automat
    {
      const T t0 = bit_notand(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x30) // function=(a and not (b)), lowered=(b notand a), set=intel
    {
      const T t0 = bit_notand(b, a);
      return t0;
    }
    else if constexpr(k ==  0x31) // function=(not (b) and (a or (c xor 1))), lowered=(b notand (a or (c xor 1))), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      const T t1 = bit_or(a, t0);
      const T t2 = bit_notand(b, t1);
      return t2;
    }
    else if constexpr(k ==  0x32) // function=(not (b) and (a or c)), lowered=(b notand (a or c)), set=automat
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_notand(b, t0);
      return t1;
    }
    else if constexpr(k ==  0x33) // function=not (b), lowered=(b xor 1), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      return t0;
    }
    else if constexpr(k ==  0x34) // function=(not ((b and c)) and (a xor b)), lowered=((b and c) notand (a xor b)), set=automat
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x35) // function=(a ? not (b) : not (c)), lowered=((b notand a) or (a notand (c xor 1))), set=intel
    {
      const T t0 = bit_notand(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_notand(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x36) // function=(b xor (a or c)), lowered=(b xor (a or c)), set=intel
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_xor(b, t0);
      return t1;
    }
    else if constexpr(k ==  0x37) // function=(b nand (a or c)), lowered=((b and (a or c)) xor 1), set=intel
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_and(b, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x38) // function=((a or c) and (a xor b)), lowered=((a or c) and (a xor b)), set=automat
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x39) // function=(b xor (a or (c xor 1))), lowered=(b xor (a or (c xor 1))), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      const T t1 = bit_or(a, t0);
      const T t2 = bit_xor(b, t1);
      return t2;
    }
    else if constexpr(k ==  0x3a) // function=(a ? not (b) : c), lowered=((b notand a) or (a notand c)), set=intel
    {
      const T t0 = bit_notand(b, a);
      const T t1 = bit_notand(a, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x3b) // function=((not (a) and c) or (b xor 1)), lowered=((a notand c) or (b xor 1)), set=automat
    {
      const T t0 = bit_notand(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x3c) // function=(b xor a), lowered=(b xor a), set=intel
    {
      const T t0 = bit_xor(b, a);
      return t0;
    }
    else if constexpr(k ==  0x3d) // function=((a xor b) or ((a or c) xor 1)), lowered=((a xor b) or ((a or c) xor 1)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T t1 = bit_or(a, c);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x3e) // function=((not (a) and c) or (a xor b)), lowered=((a notand c) or (a xor b)), set=automat
    {
      const T t0 = bit_notand(a, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x3f) // function=(b nand a), lowered=((b and a) xor 1), set=intel
    {
      const T t0 = bit_and(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0x40) // function=((not (c) and a) and b), lowered=((c notand a) and b), set=automat
    {
      const T t0 = bit_notand(c, a);
      const T t1 = bit_and(t0, b);
      return t1;
    }
    else if constexpr(k ==  0x41) // function=(c nor (b xor a)), lowered=((c or (b xor a)) xor 1), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T t1 = bit_or(c, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x42) // function=((a xor c) and (b xor c)), lowered=((a xor c) and (b xor c)), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x43) // function=(not ((a and c)) and (a xor (b xor 1))), lowered=((a and c) notand (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_and(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x44) // function=(b and not (c)), lowered=(c notand b), set=intel
    {
      const T t0 = bit_notand(c, b);
      return t0;
    }
    else if constexpr(k ==  0x45) // function=(not (c) and ((a xor 1) or b)), lowered=(c notand ((a xor 1) or b)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      const T t1 = bit_or(t0, b);
      const T t2 = bit_notand(c, t1);
      return t2;
    }
    else if constexpr(k ==  0x46) // function=(not ((a and c)) and (b xor c)), lowered=((a and c) notand (b xor c)), set=automat
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x47) // function=(b ? not (c) : not (a)), lowered=((c notand b) or (b notand (a xor 1))), set=intel
    {
      const T t0 = bit_notand(c, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_notand(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x48) // function=(b and (a xor c)), lowered=(b and (a xor c)), set=intel
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_and(b, t0);
      return t1;
    }
    else if constexpr(k ==  0x49) // function=(b ? (a xor c) : (a nor c)), lowered=((b and (a xor c)) or (b notand ((a or c) xor 1))), set=intel
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_and(b, t0);
      const T t2 = bit_or(a, c);
      const T c1 = allbits(as<T>());
      const T t3 = bit_xor(t2, c1);
      const T t4 = bit_notand(b, t3);
      const T t5 = bit_or(t1, t4);
      return t5;
    }
    else if constexpr(k ==  0x4a) // function=((b or c) and (a xor c)), lowered=((b or c) and (a xor c)), set=automat
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x4b) // function=(a xor ((b xor 1) or c)), lowered=(a xor ((b xor 1) or c)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      const T t1 = bit_or(t0, c);
      const T t2 = bit_xor(a, t1);
      return t2;
    }
    else if constexpr(k ==  0x4c) // function=(b and (a nand c)), lowered=((a and c) notand b), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_notand(t0, b);
      return t1;
    }
    else if constexpr(k ==  0x4d) // function=(b ? (a nand c) : (a nor c)), lowered=(((a and c) notand b) or (b notand ((a or c) xor 1))), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_notand(t0, b);
      const T t2 = bit_or(a, c);
      const T c1 = allbits(as<T>());
      const T t3 = bit_xor(t2, c1);
      const T t4 = bit_notand(b, t3);
      const T t5 = bit_or(t1, t4);
      return t5;
    }
    else if constexpr(k ==  0x4e) // function=(c ? not (a) : b), lowered=((a notand c) or (c notand b)), set=intel
    {
      const T t0 = bit_notand(a, c);
      const T t1 = bit_notand(c, b);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x4f) // function=((not (c) and b) or (a xor 1)), lowered=((c notand b) or (a xor 1)), set=automat
    {
      const T t0 = bit_notand(c, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x50) // function=(a and not (c)), lowered=(c notand a), set=intel
    {
      const T t0 = bit_notand(c, a);
      return t0;
    }
    else if constexpr(k ==  0x51) // function=(not (c) and (a or (b xor 1))), lowered=(c notand (a or (b xor 1))), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      const T t1 = bit_or(a, t0);
      const T t2 = bit_notand(c, t1);
      return t2;
    }
    else if constexpr(k ==  0x52) // function=(not ((b and c)) and (a xor c)), lowered=((b and c) notand (a xor c)), set=automat
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x53) // function=(a ? not (c) : not (b)), lowered=((c notand a) or (a notand (b xor 1))), set=intel
    {
      const T t0 = bit_notand(c, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_notand(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x54) // function=(not (c) and (a or b)), lowered=(c notand (a or b)), set=automat
    {
      const T t0 = bit_or(a, b);
      const T t1 = bit_notand(c, t0);
      return t1;
    }
    else if constexpr(k ==  0x55) // function=not (c), lowered=(c xor 1), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      return t0;
    }
    else if constexpr(k ==  0x56) // function=(c xor (b or a)), lowered=(c xor (b or a)), set=intel
    {
      const T t0 = bit_or(b, a);
      const T t1 = bit_xor(c, t0);
      return t1;
    }
    else if constexpr(k ==  0x57) // function=(c nand (b or a)), lowered=((c and (b or a)) xor 1), set=intel
    {
      const T t0 = bit_or(b, a);
      const T t1 = bit_and(c, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x58) // function=((a or b) and (a xor c)), lowered=((a or b) and (a xor c)), set=automat
    {
      const T t0 = bit_or(a, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x59) // function=(c xor (a or (b xor 1))), lowered=(c xor (a or (b xor 1))), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      const T t1 = bit_or(a, t0);
      const T t2 = bit_xor(c, t1);
      return t2;
    }
    else if constexpr(k ==  0x5a) // function=(c xor a), lowered=(c xor a), set=intel
    {
      const T t0 = bit_xor(c, a);
      return t0;
    }
    else if constexpr(k ==  0x5b) // function=((a xor c) or ((a or b) xor 1)), lowered=((a xor c) or ((a or b) xor 1)), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_or(a, b);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x5c) // function=(a ? not (c) : b), lowered=((c notand a) or (a notand b)), set=intel
    {
      const T t0 = bit_notand(c, a);
      const T t1 = bit_notand(a, b);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x5d) // function=((not (a) and b) or (c xor 1)), lowered=((a notand b) or (c xor 1)), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x5e) // function=((not (c) and b) or (a xor c)), lowered=((c notand b) or (a xor c)), set=automat
    {
      const T t0 = bit_notand(c, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x5f) // function=(c nand a), lowered=((c and a) xor 1), set=intel
    {
      const T t0 = bit_and(c, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0x60) // function=(a and (b xor c)), lowered=(a and (b xor c)), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_and(a, t0);
      return t1;
    }
    else if constexpr(k ==  0x61) // function=(a ? (b xor c) : (b nor c)), lowered=((a and (b xor c)) or (a notand ((b or c) xor 1))), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_and(a, t0);
      const T t2 = bit_or(b, c);
      const T c1 = allbits(as<T>());
      const T t3 = bit_xor(t2, c1);
      const T t4 = bit_notand(a, t3);
      const T t5 = bit_or(t1, t4);
      return t5;
    }
    else if constexpr(k ==  0x62) // function=((a or c) and (b xor c)), lowered=((a or c) and (b xor c)), set=automat
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x63) // function=(b xor ((a xor 1) or c)), lowered=(b xor ((a xor 1) or c)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      const T t1 = bit_or(t0, c);
      const T t2 = bit_xor(b, t1);
      return t2;
    }
    else if constexpr(k ==  0x64) // function=((a or b) and (b xor c)), lowered=((a or b) and (b xor c)), set=automat
    {
      const T t0 = bit_or(a, b);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_and(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x65) // function=(c xor ((a xor 1) or b)), lowered=(c xor ((a xor 1) or b)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      const T t1 = bit_or(t0, b);
      const T t2 = bit_xor(c, t1);
      return t2;
    }
    else if constexpr(k ==  0x66) // function=(c xor b), lowered=(c xor b), set=intel
    {
      const T t0 = bit_xor(c, b);
      return t0;
    }
    else if constexpr(k ==  0x67) // function=((b xor c) or ((a or b) xor 1)), lowered=((b xor c) or ((a or b) xor 1)), set=automat
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_or(a, b);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x68) // function=(a ? (b xor c) : (b and c)), lowered=((a and (b xor c)) or (a notand (b and c))), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_and(a, t0);
      const T t2 = bit_and(b, c);
      const T t3 = bit_notand(a, t2);
      const T t4 = bit_or(t1, t3);
      return t4;
    }
    else if constexpr(k ==  0x69) // function=(a xnor (b xor c)), lowered=((a xor (b xor c)) xor 1), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_xor(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x6a) // function=(c xor (b and a)), lowered=(c xor (b and a)), set=intel
    {
      const T t0 = bit_and(b, a);
      const T t1 = bit_xor(c, t0);
      return t1;
    }
    else if constexpr(k ==  0x6b) // function=((not (a) and c) or ((a xor 1) xor (b xor c))), lowered=((a notand c) or ((a xor 1) xor (b xor c))), set=automat
    {
      const T t0 = bit_notand(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_xor(b, c);
      const T t3 = bit_xor(t1, t2);
      const T t4 = bit_or(t0, t3);
      return t4;
    }
    else if constexpr(k ==  0x6c) // function=(b xor (a and c)), lowered=(b xor (a and c)), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_xor(b, t0);
      return t1;
    }
    else if constexpr(k ==  0x6d) // function=((not (a) and b) or ((a xor 1) xor (b xor c))), lowered=((a notand b) or ((a xor 1) xor (b xor c))), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_xor(b, c);
      const T t3 = bit_xor(t1, t2);
      const T t4 = bit_or(t0, t3);
      return t4;
    }
    else if constexpr(k ==  0x6e) // function=((not (a) and b) or (b xor c)), lowered=((a notand b) or (b xor c)), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x6f) // function=((b xor c) or (a xor 1)), lowered=((b xor c) or (a xor 1)), set=automat
    {
      const T t0 = bit_xor(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x70) // function=(a and (b nand c)), lowered=((b and c) notand a), set=intel
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_notand(t0, a);
      return t1;
    }
    else if constexpr(k ==  0x71) // function=((b nor c) or (a and (b xor c))), lowered=(((b or c) xor 1) or (a and (b xor c))), set=optimized
    {
      const T t0 = bit_or(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_xor(b, c);
      const T t3 = bit_and(a, t2);
      const T t4 = bit_or(t1, t3);
      return t4;
    }
    else if constexpr(k ==  0x72) // function=(c ? not (b) : a), lowered=((b notand c) or (c notand a)), set=intel
    {
      const T t0 = bit_notand(b, c);
      const T t1 = bit_notand(c, a);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x73) // function=((not (c) and a) or (b xor 1)), lowered=((c notand a) or (b xor 1)), set=automat
    {
      const T t0 = bit_notand(c, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x74) // function=(b ? not (c) : a), lowered=((c notand b) or (b notand a)), set=intel
    {
      const T t0 = bit_notand(c, b);
      const T t1 = bit_notand(b, a);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x75) // function=((not (b) and a) or (c xor 1)), lowered=((b notand a) or (c xor 1)), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x76) // function=((not (b) and a) or (b xor c)), lowered=((b notand a) or (b xor c)), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x77) // function=(c nand b), lowered=((c and b) xor 1), set=intel
    {
      const T t0 = bit_and(c, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0x78) // function=(a xor (b and c)), lowered=(a xor (b and c)), set=intel
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_xor(a, t0);
      return t1;
    }
    else if constexpr(k ==  0x79) // function=((not (b) and a) or ((b xor 1) xor (a xor c))), lowered=((b notand a) or ((b xor 1) xor (a xor c))), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, c);
      const T t3 = bit_xor(t1, t2);
      const T t4 = bit_or(t0, t3);
      return t4;
    }
    else if constexpr(k ==  0x7a) // function=((not (b) and a) or (a xor c)), lowered=((b notand a) or (a xor c)), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x7b) // function=((a xor c) or (b xor 1)), lowered=((a xor c) or (b xor 1)), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x7c) // function=((not (c) and a) or (a xor b)), lowered=((c notand a) or (a xor b)), set=automat
    {
      const T t0 = bit_notand(c, a);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x7d) // function=((a xor b) or (c xor 1)), lowered=((a xor b) or (c xor 1)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x7e) // function=((a xor b) or (a xor c)), lowered=((a xor b) or (a xor c)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x7f) // function=(((a and b) and c) xor 1), lowered=(((a and b) and c) xor 1), set=automat
    {
      const T t0 = bit_and(a, b);
      const T t1 = bit_and(t0, c);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x80) // function=(a and (b and c)), lowered=(a and (b and c)), set=optimized
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_and(a, t0);
      return t1;
    }
    else if constexpr(k ==  0x81) // function=(not ((a xor c)) and (a xor (b xor 1))), lowered=((a xor c) notand (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x82) // function=(c and (b xnor a)), lowered=((b xor a) notand c), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T t1 = bit_notand(t0, c);
      return t1;
    }
    else if constexpr(k ==  0x83) // function=(not ((a xor b)) and ((a xor 1) or c)), lowered=((a xor b) notand ((a xor 1) or c)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_or(t1, c);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x84) // function=(b and (a xnor c)), lowered=((a xor c) notand b), set=intel
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_notand(t0, b);
      return t1;
    }
    else if constexpr(k ==  0x85) // function=(not ((a xor c)) and (b or (c xor 1))), lowered=((a xor c) notand (b or (c xor 1))), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_or(b, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x86) // function=((b or c) and (c xor (a xor b))), lowered=((b or c) and (c xor (a xor b))), set=automat
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_xor(c, t1);
      const T t3 = bit_and(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x87) // function=(a xnor (b and c)), lowered=((a xor (b and c)) xor 1), set=intel
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_xor(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x88) // function=(c and b), lowered=(c and b), set=intel
    {
      const T t0 = bit_and(c, b);
      return t0;
    }
    else if constexpr(k ==  0x89) // function=(not ((b xor c)) and ((a xor 1) or b)), lowered=((b xor c) notand ((a xor 1) or b)), set=automat
    {
      const T t0 = bit_xor(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_or(t1, b);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x8a) // function=(not ((not (b) and a)) and c), lowered=((b notand a) notand c), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T t1 = bit_notand(t0, c);
      return t1;
    }
    else if constexpr(k ==  0x8b) // function=(b ? c : not (a)), lowered=((b and c) or (b notand (a xor 1))), set=intel
    {
      const T t0 = bit_and(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_notand(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x8c) // function=(not ((not (c) and a)) and b), lowered=((c notand a) notand b), set=automat
    {
      const T t0 = bit_notand(c, a);
      const T t1 = bit_notand(t0, b);
      return t1;
    }
    else if constexpr(k ==  0x8d) // function=(c ? b : not (a)), lowered=((c and b) or (c notand (a xor 1))), set=intel
    {
      const T t0 = bit_and(c, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_notand(c, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x8e) // function=((b and c) or (not (a) and (b xor c))), lowered=((b and c) or (a notand (b xor c))), set=optimized
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_notand(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x8f) // function=((b and c) or (a xor 1)), lowered=((b and c) or (a xor 1)), set=automat
    {
      const T t0 = bit_and(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(a, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x90) // function=(a and (b xnor c)), lowered=((b xor c) notand a), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_notand(t0, a);
      return t1;
    }
    else if constexpr(k ==  0x91) // function=(not ((b xor c)) and (a or (b xor 1))), lowered=((b xor c) notand (a or (b xor 1))), set=automat
    {
      const T t0 = bit_xor(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_or(a, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x92) // function=((a or c) and (c xor (a xor b))), lowered=((a or c) and (c xor (a xor b))), set=automat
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_xor(c, t1);
      const T t3 = bit_and(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x93) // function=(b xnor (a and c)), lowered=((b xor (a and c)) xor 1), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_xor(b, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x94) // function=((a or b) and (b xor (a xor c))), lowered=((a or b) and (b xor (a xor c))), set=automat
    {
      const T t0 = bit_or(a, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_xor(b, t1);
      const T t3 = bit_and(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x95) // function=(c xnor (b and a)), lowered=((c xor (b and a)) xor 1), set=intel
    {
      const T t0 = bit_and(b, a);
      const T t1 = bit_xor(c, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0x96) // function=(a xor (b xor c)), lowered=(a xor (b xor c)), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_xor(a, t0);
      return t1;
    }
    else if constexpr(k ==  0x97) // function=(a ? (b xnor c) : (b nand c)), lowered=(((b xor c) notand a) or (a notand ((b and c) xor 1))), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_notand(t0, a);
      const T t2 = bit_and(b, c);
      const T c1 = allbits(as<T>());
      const T t3 = bit_xor(t2, c1);
      const T t4 = bit_notand(a, t3);
      const T t5 = bit_or(t1, t4);
      return t5;
    }
    else if constexpr(k ==  0x98) // function=(not ((b xor c)) and (a or b)), lowered=((b xor c) notand (a or b)), set=automat
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_or(a, b);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0x99) // function=(c xnor b), lowered=((c xor b) xor 1), set=intel
    {
      const T t0 = bit_xor(c, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0x9a) // function=((not (b) and a) xor c), lowered=((b notand a) xor c), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T t1 = bit_xor(t0, c);
      return t1;
    }
    else if constexpr(k ==  0x9b) // function=((not (a) and c) or (b xor (c xor 1))), lowered=((a notand c) or (b xor (c xor 1))), set=automat
    {
      const T t0 = bit_notand(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x9c) // function=((not (c) and a) xor b), lowered=((c notand a) xor b), set=automat
    {
      const T t0 = bit_notand(c, a);
      const T t1 = bit_xor(t0, b);
      return t1;
    }
    else if constexpr(k ==  0x9d) // function=((not (a) and b) or (b xor (c xor 1))), lowered=((a notand b) or (b xor (c xor 1))), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x9e) // function=((b and c) or (c xor (a xor b))), lowered=((b and c) or (c xor (a xor b))), set=automat
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_xor(c, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0x9f) // function=(a nand (b xor c)), lowered=((a and (b xor c)) xor 1), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_and(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0xa0) // function=(c and a), lowered=(c and a), set=intel
    {
      const T t0 = bit_and(c, a);
      return t0;
    }
    else if constexpr(k ==  0xa1) // function=(not ((a xor c)) and (a or (b xor 1))), lowered=((a xor c) notand (a or (b xor 1))), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_or(a, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xa2) // function=(not ((not (a) and b)) and c), lowered=((a notand b) notand c), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T t1 = bit_notand(t0, c);
      return t1;
    }
    else if constexpr(k ==  0xa3) // function=(a ? c : not (b)), lowered=((a and c) or (a notand (b xor 1))), set=intel
    {
      const T t0 = bit_and(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_notand(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xa4) // function=(not ((a xor c)) and (a or b)), lowered=((a xor c) notand (a or b)), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_or(a, b);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xa5) // function=(c xnor a), lowered=((c xor a) xor 1), set=intel
    {
      const T t0 = bit_xor(c, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0xa6) // function=((not (a) and b) xor c), lowered=((a notand b) xor c), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T t1 = bit_xor(t0, c);
      return t1;
    }
    else if constexpr(k ==  0xa7) // function=((not (b) and c) or (a xor (c xor 1))), lowered=((b notand c) or (a xor (c xor 1))), set=automat
    {
      const T t0 = bit_notand(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xa8) // function=(c and (a or b)), lowered=(c and (a or b)), set=intel
    {
      const T t0 = bit_or(a, b);
      const T t1 = bit_and(c, t0);
      return t1;
    }
    else if constexpr(k ==  0xa9) // function=(c xnor (b or a)), lowered=((c xor (b or a)) xor 1), set=intel
    {
      const T t0 = bit_or(b, a);
      const T t1 = bit_xor(c, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0xaa) // function=c, lowered=c, set=intel
    {
      return c;
    }
    else if constexpr(k ==  0xab) // function=(c or (b nor a)), lowered=(c or ((b or a) xor 1)), set=intel
    {
      const T t0 = bit_or(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(c, t1);
      return t2;
    }
    else if constexpr(k ==  0xac) // function=(a ? c : b), lowered=((a and c) or (a notand b)), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_notand(a, b);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xad) // function=((b and c) or (a xor (c xor 1))), lowered=((b and c) or (a xor (c xor 1))), set=automat
    {
      const T t0 = bit_and(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xae) // function=((not (a) and b) or c), lowered=((a notand b) or c), set=automat
    {
      const T t0 = bit_notand(a, b);
      const T t1 = bit_or(t0, c);
      return t1;
    }
    else if constexpr(k ==  0xaf) // function=(c or not (a)), lowered=(c or (a xor 1)), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      const T t1 = bit_or(c, t0);
      return t1;
    }
    else if constexpr(k ==  0xb0) // function=(not ((not (c) and b)) and a), lowered=((c notand b) notand a), set=automat
    {
      const T t0 = bit_notand(c, b);
      const T t1 = bit_notand(t0, a);
      return t1;
    }
    else if constexpr(k ==  0xb1) // function=(c ? a : not (b)), lowered=((c and a) or (c notand (b xor 1))), set=intel
    {
      const T t0 = bit_and(c, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_notand(c, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xb2) // function=(b ? (a and c) : (a or c)), lowered=((b and (a and c)) or (b notand (a or c))), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_and(b, t0);
      const T t2 = bit_or(a, c);
      const T t3 = bit_notand(b, t2);
      const T t4 = bit_or(t1, t3);
      return t4;
    }
    else if constexpr(k ==  0xb3) // function=((a and c) or (b xor 1)), lowered=((a and c) or (b xor 1)), set=automat
    {
      const T t0 = bit_and(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xb4) // function=((not (c) and b) xor a), lowered=((c notand b) xor a), set=automat
    {
      const T t0 = bit_notand(c, b);
      const T t1 = bit_xor(t0, a);
      return t1;
    }
    else if constexpr(k ==  0xb5) // function=((not (b) and a) or (a xor (c xor 1))), lowered=((b notand a) or (a xor (c xor 1))), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xb6) // function=((a and c) or (c xor (a xor b))), lowered=((a and c) or (c xor (a xor b))), set=automat
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_xor(c, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xb7) // function=(b nand (a xor c)), lowered=((b and (a xor c)) xor 1), set=intel
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_and(b, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0xb8) // function=(b ? c : a), lowered=((b and c) or (b notand a)), set=intel
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_notand(b, a);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xb9) // function=((a and c) or (b xor (c xor 1))), lowered=((a and c) or (b xor (c xor 1))), set=automat
    {
      const T t0 = bit_and(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xba) // function=((not (b) and a) or c), lowered=((b notand a) or c), set=automat
    {
      const T t0 = bit_notand(b, a);
      const T t1 = bit_or(t0, c);
      return t1;
    }
    else if constexpr(k ==  0xbb) // function=(c or not (b)), lowered=(c or (b xor 1)), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      const T t1 = bit_or(c, t0);
      return t1;
    }
    else if constexpr(k ==  0xbc) // function=((a and c) or (a xor b)), lowered=((a and c) or (a xor b)), set=automat
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_xor(a, b);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xbd) // function=((a xor b) or (a xor (c xor 1))), lowered=((a xor b) or (a xor (c xor 1))), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xbe) // function=(c or (b xor a)), lowered=(c or (b xor a)), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T t1 = bit_or(c, t0);
      return t1;
    }
    else if constexpr(k ==  0xbf) // function=(c or (b nand a)), lowered=(c or ((b and a) xor 1)), set=intel
    {
      const T t0 = bit_and(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(c, t1);
      return t2;
    }
    else if constexpr(k ==  0xc0) // function=(b and a), lowered=(b and a), set=intel
    {
      const T t0 = bit_and(b, a);
      return t0;
    }
    else if constexpr(k ==  0xc1) // function=(not ((a xor b)) and (a or (c xor 1))), lowered=((a xor b) notand (a or (c xor 1))), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_or(a, t1);
      const T t3 = bit_notand(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xc2) // function=(not ((a xor b)) and (a or c)), lowered=((a xor b) notand (a or c)), set=automat
    {
      const T t0 = bit_xor(a, b);
      const T t1 = bit_or(a, c);
      const T t2 = bit_notand(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xc3) // function=(b xnor a), lowered=((b xor a) xor 1), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      return t1;
    }
    else if constexpr(k ==  0xc4) // function=(not ((not (a) and c)) and b), lowered=((a notand c) notand b), set=automat
    {
      const T t0 = bit_notand(a, c);
      const T t1 = bit_notand(t0, b);
      return t1;
    }
    else if constexpr(k ==  0xc5) // function=(a ? b : not (c)), lowered=((a and b) or (a notand (c xor 1))), set=intel
    {
      const T t0 = bit_and(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_notand(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xc6) // function=((not (a) and c) xor b), lowered=((a notand c) xor b), set=automat
    {
      const T t0 = bit_notand(a, c);
      const T t1 = bit_xor(t0, b);
      return t1;
    }
    else if constexpr(k ==  0xc7) // function=((not (c) and b) or (a xor (b xor 1))), lowered=((c notand b) or (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_notand(c, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xc8) // function=(b and (a or c)), lowered=(b and (a or c)), set=intel
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_and(b, t0);
      return t1;
    }
    else if constexpr(k ==  0xc9) // function=(b xnor (a or c)), lowered=((b xor (a or c)) xor 1), set=intel
    {
      const T t0 = bit_or(a, c);
      const T t1 = bit_xor(b, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0xca) // function=(a ? b : c), lowered=((a and b) or (a notand c)), set=intel
    {
      const T t0 = bit_and(a, b);
      const T t1 = bit_notand(a, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xcb) // function=((b and c) or (a xor (b xor 1))), lowered=((b and c) or (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_and(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xcc) // function=b, lowered=b, set=intel
    {
      return b;
    }
    else if constexpr(k ==  0xcd) // function=(b or (a nor c)), lowered=(b or ((a or c) xor 1)), set=intel
    {
      const T t0 = bit_or(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(b, t1);
      return t2;
    }
    else if constexpr(k ==  0xce) // function=((not (a) and c) or b), lowered=((a notand c) or b), set=automat
    {
      const T t0 = bit_notand(a, c);
      const T t1 = bit_or(t0, b);
      return t1;
    }
    else if constexpr(k ==  0xcf) // function=(b or not (a)), lowered=(b or (a xor 1)), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      const T t1 = bit_or(b, t0);
      return t1;
    }
    else if constexpr(k ==  0xd0) // function=(not ((not (b) and c)) and a), lowered=((b notand c) notand a), set=automat
    {
      const T t0 = bit_notand(b, c);
      const T t1 = bit_notand(t0, a);
      return t1;
    }
    else if constexpr(k ==  0xd1) // function=((b nor c) or (a and b)), lowered=(((b or c) xor 1) or (a and b)), set=optimized
    {
      const T t0 = bit_or(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_and(a, b);
      const T t3 = bit_or(t1, t2);
      return t3;
    }
    else if constexpr(k ==  0xd2) // function=((not (b) and c) xor a), lowered=((b notand c) xor a), set=automat
    {
      const T t0 = bit_notand(b, c);
      const T t1 = bit_xor(t0, a);
      return t1;
    }
    else if constexpr(k ==  0xd3) // function=((not (c) and a) or (a xor (b xor 1))), lowered=((c notand a) or (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_notand(c, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xd4) // function=((b and not (c)) or (a and (b xnor c))), lowered=((c notand b) or ((b xor c) notand a)), set=optimized
    {
      const T t0 = bit_notand(c, b);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_notand(t1, a);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xd5) // function=((a and b) or (c xor 1)), lowered=((a and b) or (c xor 1)), set=automat
    {
      const T t0 = bit_and(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xd6) // function=((a and b) or (b xor (a xor c))), lowered=((a and b) or (b xor (a xor c))), set=automat
    {
      const T t0 = bit_and(a, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_xor(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xd7) // function=(c nand (b xor a)), lowered=((c and (b xor a)) xor 1), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T t1 = bit_and(c, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0xd8) // function=(c ? b : a), lowered=((c and b) or (c notand a)), set=intel
    {
      const T t0 = bit_and(c, b);
      const T t1 = bit_notand(c, a);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xd9) // function=((a and b) or (b xor (c xor 1))), lowered=((a and b) or (b xor (c xor 1))), set=automat
    {
      const T t0 = bit_and(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(b, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xda) // function=((a and b) or (a xor c)), lowered=((a and b) or (a xor c)), set=automat
    {
      const T t0 = bit_and(a, b);
      const T t1 = bit_xor(a, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xdb) // function=((a xor c) or (a xor (b xor 1))), lowered=((a xor c) or (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_xor(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xdc) // function=((not (c) and a) or b), lowered=((c notand a) or b), set=automat
    {
      const T t0 = bit_notand(c, a);
      const T t1 = bit_or(t0, b);
      return t1;
    }
    else if constexpr(k ==  0xdd) // function=(b or not (c)), lowered=(b or (c xor 1)), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      const T t1 = bit_or(b, t0);
      return t1;
    }
    else if constexpr(k ==  0xde) // function=(b or (a xor c)), lowered=(b or (a xor c)), set=intel
    {
      const T t0 = bit_xor(a, c);
      const T t1 = bit_or(b, t0);
      return t1;
    }
    else if constexpr(k ==  0xdf) // function=(b or (a nand c)), lowered=(b or ((a and c) xor 1)), set=intel
    {
      const T t0 = bit_and(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(b, t1);
      return t2;
    }
    else if constexpr(k ==  0xe0) // function=(a and (b or c)), lowered=(a and (b or c)), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_and(a, t0);
      return t1;
    }
    else if constexpr(k ==  0xe1) // function=(a xnor (b or c)), lowered=((a xor (b or c)) xor 1), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_xor(a, t0);
      const T c1 = allbits(as<T>());
      const T t2 = bit_xor(t1, c1);
      return t2;
    }
    else if constexpr(k ==  0xe2) // function=(b ? a : c), lowered=((b and a) or (b notand c)), set=intel
    {
      const T t0 = bit_and(b, a);
      const T t1 = bit_notand(b, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xe3) // function=((a and c) or (a xor (b xor 1))), lowered=((a and c) or (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_and(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xe4) // function=(c ? a : b), lowered=((c and a) or (c notand b)), set=intel
    {
      const T t0 = bit_and(c, a);
      const T t1 = bit_notand(c, b);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xe5) // function=((a and b) or (a xor (c xor 1))), lowered=((a and b) or (a xor (c xor 1))), set=automat
    {
      const T t0 = bit_and(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xe6) // function=((a and b) or (b xor c)), lowered=((a and b) or (b xor c)), set=automat
    {
      const T t0 = bit_and(a, b);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_or(t0, t1);
      return t2;
    }
    else if constexpr(k ==  0xe7) // function=((b xor c) or (a xor (b xor 1))), lowered=((b xor c) or (a xor (b xor 1))), set=automat
    {
      const T t0 = bit_xor(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(b, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xe8) // function=((b and c) or (a and (b xor c))), lowered=((b and c) or (a and (b xor c))), set=optimized
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_xor(b, c);
      const T t2 = bit_and(a, t1);
      const T t3 = bit_or(t0, t2);
      return t3;
    }
    else if constexpr(k ==  0xe9) // function=((a and b) or (b xor (a xor (c xor 1)))), lowered=((a and b) or (b xor (a xor (c xor 1)))), set=automat
    {
      const T t0 = bit_and(a, b);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(c, c1);
      const T t2 = bit_xor(a, t1);
      const T t3 = bit_xor(b, t2);
      const T t4 = bit_or(t0, t3);
      return t4;
    }
    else if constexpr(k ==  0xea) // function=(c or (b and a)), lowered=(c or (b and a)), set=intel
    {
      const T t0 = bit_and(b, a);
      const T t1 = bit_or(c, t0);
      return t1;
    }
    else if constexpr(k ==  0xeb) // function=(c or (b xnor a)), lowered=(c or ((b xor a) xor 1)), set=intel
    {
      const T t0 = bit_xor(b, a);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(c, t1);
      return t2;
    }
    else if constexpr(k ==  0xec) // function=(b or (a and c)), lowered=(b or (a and c)), set=intel
    {
      const T t0 = bit_and(a, c);
      const T t1 = bit_or(b, t0);
      return t1;
    }
    else if constexpr(k ==  0xed) // function=(b or (a xnor c)), lowered=(b or ((a xor c) xor 1)), set=intel
    {
      const T t0 = bit_xor(a, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(b, t1);
      return t2;
    }
    else if constexpr(k ==  0xee) // function=(c or b), lowered=(c or b), set=intel
    {
      const T t0 = bit_or(c, b);
      return t0;
    }
    else if constexpr(k ==  0xef) // function=(b or ((a xor 1) or c)), lowered=(b or ((a xor 1) or c)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(a, c1);
      const T t1 = bit_or(t0, c);
      const T t2 = bit_or(b, t1);
      return t2;
    }
    else if constexpr(k ==  0xf0) // function= a, lowered= a, set=intel
    {
      return a;
    }
    else if constexpr(k ==  0xf1) // function=(a or (b nor c)), lowered=(a or ((b or c) xor 1)), set=intel
    {
      const T t0 = bit_or(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(a, t1);
      return t2;
    }
    else if constexpr(k ==  0xf2) // function=((not (b) and c) or a), lowered=((b notand c) or a), set=automat
    {
      const T t0 = bit_notand(b, c);
      const T t1 = bit_or(t0, a);
      return t1;
    }
    else if constexpr(k ==  0xf3) // function=(a or not (b)), lowered=(a or (b xor 1)), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      const T t1 = bit_or(a, t0);
      return t1;
    }
    else if constexpr(k ==  0xf4) // function=((not (c) and b) or a), lowered=((c notand b) or a), set=automat
    {
      const T t0 = bit_notand(c, b);
      const T t1 = bit_or(t0, a);
      return t1;
    }
    else if constexpr(k ==  0xf5) // function=(a or not (c)), lowered=(a or (c xor 1)), set=intel
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      const T t1 = bit_or(a, t0);
      return t1;
    }
    else if constexpr(k ==  0xf6) // function=(a or (b xor c)), lowered=(a or (b xor c)), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T t1 = bit_or(a, t0);
      return t1;
    }
    else if constexpr(k ==  0xf7) // function=(a or (b nand c)), lowered=(a or ((b and c) xor 1)), set=intel
    {
      const T t0 = bit_and(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(a, t1);
      return t2;
    }
    else if constexpr(k ==  0xf8) // function=(a or (b and c)), lowered=(a or (b and c)), set=intel
    {
      const T t0 = bit_and(b, c);
      const T t1 = bit_or(a, t0);
      return t1;
    }
    else if constexpr(k ==  0xf9) // function=(a or (b xnor c)), lowered=(a or ((b xor c) xor 1)), set=intel
    {
      const T t0 = bit_xor(b, c);
      const T c1 = allbits(as<T>());
      const T t1 = bit_xor(t0, c1);
      const T t2 = bit_or(a, t1);
      return t2;
    }
    else if constexpr(k ==  0xfa) // function=(c or a), lowered=(c or a), set=intel
    {
      const T t0 = bit_or(c, a);
      return t0;
    }
    else if constexpr(k ==  0xfb) // function=(a or ((b xor 1) or c)), lowered=(a or ((b xor 1) or c)), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(b, c1);
      const T t1 = bit_or(t0, c);
      const T t2 = bit_or(a, t1);
      return t2;
    }
    else if constexpr(k ==  0xfc) // function=(b or a), lowered=(b or a), set=intel
    {
      const T t0 = bit_or(b, a);
      return t0;
    }
    else if constexpr(k ==  0xfd) // function=(a or (b or (c xor 1))), lowered=(a or (b or (c xor 1))), set=automat
    {
      const T c1 = allbits(as<T>());
      const T t0 = bit_xor(c, c1);
      const T t1 = bit_or(b, t0);
      const T t2 = bit_or(a, t1);
      return t2;
    }
    else if constexpr(k ==  0xfe) // function=(a or (b or c)), lowered=(a or (b or c)), set=intel
    {
      const T t0 = bit_or(b, c);
      const T t1 = bit_or(a, t0);
      return t1;
    }
    else if constexpr(k ==  0xff) // function=1, lowered=1, set=intel
    {
      const T c1 = allbits(as<T>());
      return c1;
    }
  }

}
