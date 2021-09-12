//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include "test.hpp"
#include <eve/concept/value.hpp>
#include <eve/constant/valmin.hpp>
#include <eve/constant/valmax.hpp>
#include <eve/function/all.hpp>
#include <eve/function/cyl_bessel_jn.hpp>
#include <eve/function/is_negative.hpp>
#include <eve/function/is_positive.hpp>
#include <type_traits>
#include <cmath>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/constant/smallestposval.hpp>
#include <eve/platform.hpp>
#include <cmath>

//==================================================================================================
// Types tests
//==================================================================================================
// EVE_TEST_TYPES( "Check return types of cyl_bessel_jn"
//             , eve::test::simd::ieee_reals
//             )
// <typename T>(eve::as<T>)
// {
//   using v_t = eve::element_type_t<T>;
//   using i_t = eve::as_integer_t<v_t>;
//   using I_t = eve::as_integer_t<T>;

//   TTS_EXPR_IS( eve::cyl_bessel_jn(T(), T())  ,   T);
//   TTS_EXPR_IS( eve::cyl_bessel_jn(v_t(),v_t()), v_t);
//   TTS_EXPR_IS( eve::cyl_bessel_jn(i_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_jn(I_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_jn(I_t(),v_t()), T);
//   TTS_EXPR_IS( eve::cyl_bessel_jn(i_t(),v_t()), v_t);
// };

//==================================================================================================
// cyl_bessel_jn  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_jn on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0),
                              eve::test::ramp(0))
        )
<typename T>(T const& a0, T const & n)
{
  using v_t = eve::element_type_t<T>;
  auto eve__cyl_bessel_jn =  [](auto n, auto x) {
    auto j = eve::cyl_bessel_jn(n, x) ;
    return kumi::get<0>(j);
  };
  auto boost_cyl_bessel_jn =  [](v_t n, v_t x) { return std::cyl_bessel_j(n, x); };
  if constexpr( eve::platform::supports_invalids )
  {
//     TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(2), eve::minf(eve::as<v_t>())), eve::zero(eve::as<v_t>()), 0);
//     TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(2), eve::inf(eve::as<v_t>())), eve::zero(eve::as<v_t>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(3), eve::nan(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
  }
//   if constexpr(eve::cardinal_v<T> == 1)
//   {
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(0), v_t(0)), v_t(boost_cyl_bessel_jn(0.0, 0.0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(1), v_t(0)), v_t(boost_cyl_bessel_jn(1.0, 0.0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(2), v_t(0)), v_t(boost_cyl_bessel_jn(2.0, 0.0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(0), v_t(1)), v_t(boost_cyl_bessel_jn(0.0, 1.0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(1), v_t(1)), v_t(boost_cyl_bessel_jn(1.0, 1.0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(2), v_t(1)), v_t(boost_cyl_bessel_jn(2.0, 1.0)), 10);
//     for(int i=1; i < 2; i*= 2)
//     {
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(i), v_t(10)), v_t(boost_cyl_bessel_jn(i, 10)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(i), v_t(5)), v_t(boost_cyl_bessel_jn(i, 5)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(i), v_t(2)), v_t(boost_cyl_bessel_jn(i, 2)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(i), v_t(1)), v_t(boost_cyl_bessel_jn(i, 1)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(v_t(i), v_t(0)), v_t(boost_cyl_bessel_jn(i, 0)), 10);
//     }
//   }
//      TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(0), T(0)), T(boost_cyl_bessel_jn(0.0, 0.0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(1), T(0)), T(boost_cyl_bessel_jn(1.0, 0.0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(2), T(0)), T(boost_cyl_bessel_jn(2.0, 0.0)), 10);
//        TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(0), T(1)), T(boost_cyl_bessel_jn(0.0, 1.0)), 10);
//        TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(1), T(1)), T(boost_cyl_bessel_jn(1.0, 1.0)), 10);
//        TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(2), T(1)), T(boost_cyl_bessel_jn(2.0, 1.0)), 10);
//        TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(0), T(3)), T(boost_cyl_bessel_jn(0.0, 3.0)), 10);
//        TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(1), T(3)), T(boost_cyl_bessel_jn(1.0, 3.0)), 10);
//        TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(2), T(3)), T(boost_cyl_bessel_jn(2.0, 3.0)), 10);

//   TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(0), a0), map(boost_cyl_bessel_jn, T(0), a0), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(1), a0), map(boost_cyl_bessel_jn, T(1), a0), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(2), a0), map(boost_cyl_bessel_jn, T(2), a0), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_jn(T(0.5), a0), map(boost_cyl_bessel_jn, T(0.5), a0), 10);
  TTS_ULP_EQUAL(eve__cyl_bessel_jn(n     , a0), map(boost_cyl_bessel_jn, n,    a0), 100);

};
