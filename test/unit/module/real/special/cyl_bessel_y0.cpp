//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include <eve/function/cyl_bessel_y0.hpp>
#include <eve/function/prev.hpp>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/platform.hpp>
#include <cmath>

EVE_TEST_TYPES( "Check return types of cyl_bessel_y0"
            , eve::test::simd::ieee_reals
            )
<typename T>(eve::as<T>)
{
  using v_t = eve::element_type_t<T>;
  TTS_EXPR_IS(eve::cyl_bessel_y0(T(0)), T);
  TTS_EXPR_IS(eve::cyl_bessel_y0(v_t(0)), v_t);
};

 EVE_TEST( "Check behavior of cyl_bessel_y0 on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0))
         )
   <typename T>(T const& a0)
{
  using v_t = eve::element_type_t<T>;

  auto eve__cyl_bessel_y0 =  [](auto x) { return eve::cyl_bessel_y0(x); };
  auto std__cyl_bessel_y0 =  [](auto x)->v_t { return std::cyl_neumann(0, x); };
  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(eve__cyl_bessel_y0(eve::minf(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y0(eve::inf(eve::as<v_t>())), v_t(0), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y0(eve::nan(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y0(eve::minf(eve::as< T>())), eve::nan(eve::as< T>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y0(eve::inf(eve::as< T>())),  T(0), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_y0(eve::nan(eve::as< T>())), eve::nan(eve::as< T>()), 0);
  }
  TTS_ULP_EQUAL(eve__cyl_bessel_y0(v_t(10)), v_t(5.567116728359938e-02), 1.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0(v_t(5)), v_t(-3.085176252490338e-01), 1.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0(v_t(2)), v_t(  5.103756726497453e-01), 1.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0(v_t(1.5)),v_t(  3.824489237977589e-01 ), 2.0);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0(v_t(0.5)),v_t( -4.445187335067066e-01 ), 2.0);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0(v_t(1)), v_t( 8.825696421567700e-02), 2.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0(v_t(0)), eve::minf(eve::as<v_t>()), 0.0);

  TTS_ULP_EQUAL(eve__cyl_bessel_y0( T(10)),  T(5.567116728359938e-02), 1.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0( T(5)),  T(-3.085176252490338e-01), 1.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0( T(2)),  T(  5.103756726497453e-01), 1.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0( T(1.5)), T(  3.824489237977589e-01 ), 2.0);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0( T(0.5)), T( -4.445187335067066e-01 ), 2.0);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0( T(1)),  T( 8.825696421567700e-02), 2.5);
  TTS_ULP_EQUAL(eve__cyl_bessel_y0( T(0)), eve::minf(eve::as< T>()), 0.0);


  TTS_ULP_EQUAL(eve__cyl_bessel_y0(a0), map(std__cyl_bessel_y0, a0), 100.0);
};
