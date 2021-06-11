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
#include <eve/function/cyl_bessel_j0.hpp>
#include <eve/function/diff/cyl_bessel_j0.hpp>
#include <eve/function/is_negative.hpp>
#include <eve/function/is_positive.hpp>
#include <type_traits>
#include <cmath>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/constant/smallestposval.hpp>
#include <eve/platform.hpp>

//==================================================================================================
// Types tests
//==================================================================================================
EVE_TEST_TYPES( "Check return types of cyl_bessel_j0"
            , eve::test::simd::ieee_reals
            )
<typename T>(eve::as_<T>)
{
  using v_t = eve::element_type_t<T>;

  TTS_EXPR_IS( eve::cyl_bessel_j0(T())  , T);
  TTS_EXPR_IS( eve::cyl_bessel_j0(v_t()), v_t);
};

//==================================================================================================
// cyl_bessel_j0  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_j0 on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0))
        )
<typename T>(T const& a0 )
{
  using v_t = eve::element_type_t<T>;
  using eve::cyl_bessel_j0;
  using eve::as;
  TTS_ULP_EQUAL( cyl_bessel_j0(a0),  map([&](auto e) -> v_t{ return std::cyl_bessel_j(0, e); }, a0), 23);
  auto dcyl_bessel_j0 = [](auto e) -> v_t{return  -std::cyl_bessel_j(1, e);};
  TTS_ULP_EQUAL( eve::diff(cyl_bessel_j0)(a0),  map(dcyl_bessel_j0, a0), 4);

  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(cyl_bessel_j0(eve::minf(eve::as<T>())), T(0), 0);
    TTS_ULP_EQUAL(cyl_bessel_j0(eve::inf(eve::as<T>())), T(0), 0);
    TTS_ULP_EQUAL(cyl_bessel_j0(eve::nan(eve::as<T>())), eve::nan(eve::as<T>()), 0);
  }
  TTS_ULP_EQUAL(cyl_bessel_j0(T(10)), T( -2.459357644513482e-01), 2.5);
  TTS_ULP_EQUAL(cyl_bessel_j0(T(5)), T(-1.775967713143384e-01), 1.5);
  TTS_ULP_EQUAL(cyl_bessel_j0(T(2)), T( 2.238907791412356e-01), 0.5);
  TTS_ULP_EQUAL(cyl_bessel_j0(T(1)), T( 7.651976865579666e-01), 0.5);
  TTS_ULP_EQUAL(cyl_bessel_j0(T(0)), T(1), 0);
  TTS_ULP_EQUAL(cyl_bessel_j0(T(0.5)), T(std::cyl_bessel_j(0, v_t(0.5))), 1. );
};
