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
#include <eve/function/cyl_bessel_j1.hpp>
#include <type_traits>
#include <cmath>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/platform.hpp>

//==================================================================================================
// Types tests
//==================================================================================================
EVE_TEST_TYPES( "Check return types of cyl_bessel_j1"
            , eve::test::simd::ieee_reals
            )
<typename T>(eve::as_<T>)
{
  using v_t = eve::element_type_t<T>;

  TTS_EXPR_IS( eve::cyl_bessel_j1(T())  , T);
  TTS_EXPR_IS( eve::cyl_bessel_j1(v_t()), v_t);
};

//==================================================================================================
// cyl_bessel_j1  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_j1 on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0))
        )
<typename T>(T const& a0 )
{
  using v_t = eve::element_type_t<T>;
  using eve::cyl_bessel_j1;
  using eve::as;
  TTS_ULP_EQUAL( cyl_bessel_j1(a0),  map([&](auto e) -> v_t{ return std::cyl_bessel_j(1, e); }, a0), 4);

  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(cyl_bessel_j1(eve::minf(eve::as<T>())), T(0), 0);
    TTS_ULP_EQUAL(cyl_bessel_j1(eve::inf(eve::as<T>())), T(0), 0);
    TTS_ULP_EQUAL(cyl_bessel_j1(eve::nan(eve::as<T>())), eve::nan(eve::as<T>()), 0);
  }
  TTS_ULP_EQUAL(cyl_bessel_j1(T(20)), T( 6.683312417584991e-02), 11);
  TTS_ULP_EQUAL(cyl_bessel_j1(T(10)), T( 4.347274616886149e-02), 8);
  TTS_ULP_EQUAL(cyl_bessel_j1(T(5)), T( -3.275791375914651e-01), 1.5);
  TTS_ULP_EQUAL(cyl_bessel_j1(T(2)), T( 5.767248077568736e-01), 1);
  TTS_ULP_EQUAL(cyl_bessel_j1(T(1)), T(4.400505857449336e-01), 0.5);
  TTS_ULP_EQUAL(cyl_bessel_j1(T(0)), T(0), 0);
  TTS_ULP_EQUAL(cyl_bessel_j1(T(0.5)), T(std::cyl_bessel_j(1, v_t(0.5))), 1. );
};
