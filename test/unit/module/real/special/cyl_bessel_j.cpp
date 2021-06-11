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
#include <eve/function/cyl_bessel_j.hpp>
#include <eve/function/diff/cyl_bessel_j.hpp>
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
EVE_TEST_TYPES( "Check return types of cyl_bessel_j"
            , eve::test::simd::ieee_reals
            )
<typename T>(eve::as_<T>)
{
  using v_t = eve::element_type_t<T>;
  using i_t = eve::as_integer_t<v_t>;
//  using I_t = eve::as_integer_t<T>;

  TTS_EXPR_IS( eve::cyl_bessel_j(T(), T())  ,   T);
  TTS_EXPR_IS( eve::cyl_bessel_j(v_t(),v_t()), v_t);
  TTS_EXPR_IS( eve::cyl_bessel_j(i_t(),T()),   T);
  // TTS_EXPR_IS( eve::cyl_bessel_j(I_t(),T()),   T);
  //  TTS_EXPR_IS( eve::cyl_bessel_j(I_t(),v_t()), T);
  TTS_EXPR_IS( eve::cyl_bessel_j(i_t(),v_t()), v_t);
};

//==================================================================================================
// cyl_bessel_j  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_j on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0))
        )
<typename T>(T const& a0 )
{
  using v_t = eve::element_type_t<T>;
  using eve::cyl_bessel_j;
  using eve::as;
  for(int i=1; i < 4 ; ++i)
  {
    TTS_ULP_EQUAL( cyl_bessel_j(i, a0),  map([i](auto e) -> v_t {return std::cyl_bessel_j(i, e);}, a0), 100);
    auto dcyl_bessel_j = [i](auto e)-> v_t {return -cyl_bessel_j(i+1, e)+cyl_bessel_j(i, e)*i/e;};
    TTS_ULP_EQUAL( eve::diff(cyl_bessel_j)(i, a0),  map(dcyl_bessel_j, a0), 25);
  }

  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(cyl_bessel_j(T(2), eve::minf(eve::as<T>())), eve::zero(eve::as<T>()), 0);
    TTS_ULP_EQUAL(cyl_bessel_j(T(2), eve::inf(eve::as<T>())), eve::zero(eve::as<T>()), 0);
    TTS_ULP_EQUAL(cyl_bessel_j(T(3), eve::nan(eve::as<T>())), eve::nan(eve::as<T>()), 0);
  }

  for(int i=1; i < 5; ++i)
  {
    TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(10)), T(std::cyl_bessel_j(i, 10)), 20);
    TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(5)), T(std::cyl_bessel_j(i, 5)), 20);
    TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(2)), T(std::cyl_bessel_j(i, 2)), 20);
    TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(1)), T(std::cyl_bessel_j(i, 1)), 20);
    TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(0)), T(std::cyl_bessel_j(i, 0)), 20);
    TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(0.5)), T(std::cyl_bessel_j(i, 0.5)), 20);
  }

  if constexpr( eve::platform::supports_invalids )
  {
    TTS_IEEE_EQUAL(cyl_bessel_j(T(1), eve::nan(eve::as<T>()))  , eve::nan(eve::as<T>()) );
    TTS_IEEE_EQUAL(cyl_bessel_j(T(1), eve::inf(eve::as<T>()))   , T(0) );
  }

};
