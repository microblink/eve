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
#include <eve/function/cyl_bessel_y0.hpp>
#include <type_traits>
#include <cmath>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/platform.hpp>

//==================================================================================================
// Types tests
//==================================================================================================
EVE_TEST_TYPES( "Check return types of cyl_bessel_y0"
            , eve::test::simd::ieee_reals
            )
<typename T>(eve::as_<T>)
{
  using v_t = eve::element_type_t<T>;

  TTS_EXPR_IS( eve::cyl_bessel_y0(T())  , T);
  TTS_EXPR_IS( eve::cyl_bessel_y0(v_t()), v_t);
};

//==================================================================================================
// cyl_bessel_y0  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_y0 on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(1.0, 4.0))
        )
<typename T>(T const& a0 )
{
  using v_t = eve::element_type_t<T>;
  using eve::cyl_bessel_y0;
  using eve::as;
  TTS_ULP_EQUAL( cyl_bessel_y0(a0),  map([&](auto e) -> v_t{ return std::cyl_neumann(0, v_t(e)); }, a0), 4);

  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(cyl_bessel_y0(eve::minf(eve::as<T>())), eve::nan(eve::as<T>()), 0);
    TTS_ULP_EQUAL(cyl_bessel_y0(eve::inf(eve::as<T>())), T(0), 0);
    TTS_ULP_EQUAL(cyl_bessel_y0(eve::nan(eve::as<T>())), eve::nan(eve::as<T>()), 0);
  }
  TTS_ULP_EQUAL(cyl_bessel_y0(T(20)), T(std::cyl_neumann(0, v_t(20))), 20.0);
  TTS_ULP_EQUAL(cyl_bessel_y0(T(10)), T(std::cyl_neumann(0, v_t(10))), 10.0);
  TTS_ULP_EQUAL(cyl_bessel_y0(T(5)) , T(std::cyl_neumann(0, v_t(5))) , 2.0);
  TTS_ULP_EQUAL(cyl_bessel_y0(T(2)) , T(std::cyl_neumann(0, v_t(2))) , 120.0);
  TTS_ULP_EQUAL(cyl_bessel_y0(T(1)) , T(std::cyl_neumann(0, v_t(1))) , 2.0);
  TTS_ULP_EQUAL(cyl_bessel_y0(T(0)) , T(std::cyl_neumann(0, v_t(0))) , 2.0);
  TTS_ULP_EQUAL(cyl_bessel_y0(T(0.5)), T(std::cyl_neumann(0, v_t(0.5))), 2.9 );
};
