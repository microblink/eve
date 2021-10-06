// //==================================================================================================
// /**
//   EVE - Expressive Vector Engine
//   Copyright : EVE Contributors & Maintainers
//   SPDX-License-Identifier: MIT
// **/
// //==================================================================================================
#include "test.hpp"
#include <eve/concept/value.hpp>
#include <eve/constant/valmin.hpp>
#include <eve/constant/valmax.hpp>
#include <eve/function/all.hpp>
#include <eve/function/cyl_bessel_yn.hpp>
#include <eve/function/cyl_bessel_y1.hpp>
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
#include <boost/math/special_functions/bessel.hpp>

//==================================================================================================
// Types tests
//==================================================================================================
// EVE_TEST_TYPES( "Check return types of cyl_bessel_jy"
//             , eve::test::simd::ieee_reals
//             )
// <typename T>(eve::as<T>)
// {
//   using v_t = eve::element_type_t<T>;
//   using i_t = eve::as_integer_t<v_t>;
//   using I_t = eve::as_integer_t<T>;
//   TTS_EXPR_IS( eve::cyl_bessel_yn(T(), T())  ,   T);
//   TTS_EXPR_IS( eve::cyl_bessel_yn(v_t(),v_t()), v_t);
//   TTS_EXPR_IS( eve::cyl_bessel_yn(i_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_yn(I_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_yn(I_t(),v_t()), T);
//   TTS_EXPR_IS( eve::cyl_bessel_yn(i_t(),v_t()), v_t);
// };

//==================================================================================================
//===  cyl_bessel_yn scalar tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_yn for j on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0),
                              eve::test::ramp(0))
        )
  <typename T>(T const& a0, T const & n)
{
  using v_t = eve::element_type_t<T>;
  auto eve__cyl_bessel_yn =  [](auto n, auto x) {
    return eve::cyl_bessel_yn(n, x) ;
  };
  auto std_cyl_bessel_yn =  [](v_t n, v_t x) ->v_t{ return /*boost::math*/std::cyl_neumann(n, x); };
  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(2), eve::minf(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(2), eve::inf(eve::as<v_t>())), eve::zero(eve::as<v_t>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(3), eve::nan(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
  }
  if constexpr(eve::cardinal_v<T> == 1)
  {
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(0), v_t(0)), v_t(std_cyl_bessel_yn(v_t(0), v_t(0))), 10);
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(1), v_t(0)), v_t(std_cyl_bessel_yn(v_t(1), v_t(0))), 10);
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(2), v_t(0)), v_t(std_cyl_bessel_yn(v_t(2), v_t(0))), 10);
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(0), v_t(1)), v_t(std_cyl_bessel_yn(v_t(0), v_t(1))), 10);
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(1), v_t(1)), v_t(std_cyl_bessel_yn(v_t(1), v_t(1))), 10);
    TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(2), v_t(1)), v_t(std_cyl_bessel_yn(v_t(2), v_t(1))), 10);
    for(int i=3; i < 4; ++i)
    {
      TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(i), v_t(10)), v_t(std_cyl_bessel_yn(v_t(i), v_t(10))), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(i), v_t(5)), v_t(std_cyl_bessel_yn(v_t(i), v_t(5))), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(i), v_t(2)), v_t(std_cyl_bessel_yn(v_t(i), v_t(2))), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(i), v_t(1)), v_t(std_cyl_bessel_yn(v_t(i), v_t(1))), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_yn(v_t(i), v_t(0)), v_t(std_cyl_bessel_yn(v_t(i), v_t(0))), 10);
    }
  }

  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(0), T(0)), T(std_cyl_bessel_yn(v_t(0.0), v_t(0.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(1), T(0)), T(std_cyl_bessel_yn(v_t(1.0), v_t(0.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(2), T(0)), T(std_cyl_bessel_yn(v_t(2.0), v_t(0.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(0), T(1)), T(std_cyl_bessel_yn(v_t(0.0), v_t(1.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(1), T(1)), T(std_cyl_bessel_yn(v_t(1.0), v_t(1.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(2), T(1)), T(std_cyl_bessel_yn(v_t(2.0), v_t(1.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(0), T(3)), T(std_cyl_bessel_yn(v_t(0.0), v_t(3.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(1), T(3)), T(std_cyl_bessel_yn(v_t(1.0), v_t(3.0))), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(2), T(3)), T(std_cyl_bessel_yn(v_t(2.0), v_t(3.0))), 200);

  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(0), a0), map(std_cyl_bessel_yn, v_t(0), a0), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(1), a0), map(std_cyl_bessel_yn, v_t(1), a0), 200);
  TTS_ULP_EQUAL(eve::cyl_bessel_y1(a0),       map(std_cyl_bessel_yn, v_t(1), a0), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(T(2), a0), map(std_cyl_bessel_yn, v_t(2), a0), 200);
  TTS_ULP_EQUAL(eve__cyl_bessel_yn(n   , a0), map(std_cyl_bessel_yn, n     , a0), 200);
};
