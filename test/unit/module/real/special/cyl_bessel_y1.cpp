//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include <eve/function/cyl_bessel_y1.hpp>
#include <eve/constant/eps.hpp>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/platform.hpp>
#include <cmath>

EVE_TEST_TYPES( "Check return types of cyl_bessel_y1"
              , eve::test::simd::ieee_reals
              )
  <typename T>(eve::as<T>)
{
  using v_t = eve::element_type_t<T>;
  TTS_EXPR_IS(eve::cyl_bessel_y1(T(0)), T);
  TTS_EXPR_IS(eve::cyl_bessel_y1(v_t(0)), v_t);
};

EVE_TEST( "Check behavior of cyl_bessel_y1 on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0))
        )
  <typename T>(T const& a0)
{
  using v_t = eve::element_type_t<T>;

  auto eve__cyl_bessel_y1 =  [](auto x) { return eve::cyl_bessel_y1(x); };
  auto std__cyl_bessel_y1 =  [](auto x)->v_t { return std::cyl_neumann(1, x); };
  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(eve__cyl_bessel_y1(eve::minf(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y1(eve::inf(eve::as<v_t>())), v_t(0), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y1(eve::nan(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);

    TTS_ULP_EQUAL(eve__cyl_bessel_y1(eve::minf(eve::as<T>())), eve::nan(eve::as<T>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y1(eve::inf(eve::as<T>())), T(0), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y1(eve::nan(eve::as<T>())), eve::nan(eve::as<T>()), 0);
  }

  TTS_ULP_EQUAL(eve__cyl_bessel_y1(v_t(10)), v_t(  2.490154242069539e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(v_t(5)), v_t(1.478631433912269e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(v_t(2)), v_t( -1.070324315409375e-01), 4);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(v_t(1.5)), v_t(  -4.123086269739113e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(v_t(1)), v_t( -7.812128213002889e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(v_t(0)), eve::minf(eve::as<v_t>()), 0);

  TTS_ULP_EQUAL(eve__cyl_bessel_y1(T(10)), T(  2.490154242069539e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(T(5)), T(1.478631433912269e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(T(2)), T( -1.070324315409375e-01), 4);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(T(1.5)), T(  -4.123086269739113e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(T(1)), T( -7.812128213002889e-01), 2);
  TTS_ULP_EQUAL(eve__cyl_bessel_y1(T(0)), eve::minf(eve::as<T>()), 0);


  TTS_RELATIVE_EQUAL(eve__cyl_bessel_y1(a0), map(std__cyl_bessel_y1, a0), 0.001);
};
