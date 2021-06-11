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
#include <eve/function/cyl_bessel_y.hpp>
#include <eve/function/diff/cyl_bessel_y.hpp>
#include <type_traits>
#include <cmath>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/constant/smallestposval.hpp>
#include <eve/platform.hpp>
#include <boost/math/special_functions/bessel.hpp>

//==================================================================================================
// Types tests
//==================================================================================================
EVE_TEST_TYPES( "Check return types of cyl_bessel_y"
            , eve::test::simd::ieee_reals
            )
<typename T>(eve::as_<T>)
{
  using v_t = eve::element_type_t<T>;
  using i_t = eve::as_integer_t<v_t>;
//  using I_t = eve::as_integer_t<T>;

  TTS_EXPR_IS( eve::cyl_bessel_y(T(), T())  ,   T);
  TTS_EXPR_IS( eve::cyl_bessel_y(v_t(),v_t()), v_t);
  TTS_EXPR_IS( eve::cyl_bessel_y(i_t(),T()),   T);
  // TTS_EXPR_IS( eve::cyl_bessel_y(I_t(),T()),   T);
  //  TTS_EXPR_IS( eve::cyl_bessel_y(I_t(),v_t()), T);
  TTS_EXPR_IS( eve::cyl_bessel_y(i_t(),v_t()), v_t);
};

//==================================================================================================
// cyl_bessel_y  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_y on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 5.0))
        )
<typename T>(T const& a0 )
{
  using v_t = eve::element_type_t<T>;
  using eve::cyl_bessel_y;
  using eve::as;
  for(int i=1; i < 4 ; ++i)
  {
    std::cout << "i =  "<< i << std::endl;
    TTS_ULP_EQUAL( cyl_bessel_y(i, a0),  map([i](auto e) -> v_t {return boost::math::cyl_neumann(i, v_t(e));}, a0), 4000);

    for(int j=0; j < eve::cardinal_v<T>; ++j)
    {
      std::cout << "j =  " << j << std::endl;
      TTS_ULP_EQUAL( cyl_bessel_y(i, a0.get(j)),  boost::math::cyl_neumann(i, a0.get(j)), 4000);
    }
    auto dcyl_bessel_y = [i](auto e)-> v_t {return -cyl_bessel_y(i+1, e)+cyl_bessel_y(i, e)*i/e;};
    TTS_ULP_EQUAL( eve::diff(cyl_bessel_y)(i, a0),  map(dcyl_bessel_y, a0), 25);
  }

  if constexpr( eve::platform::supports_invalids )
  {
    TTS_ULP_EQUAL(cyl_bessel_y(T(2), eve::inf(eve::as<T>())), eve::zero(eve::as<T>()), 0);
    TTS_ULP_EQUAL(cyl_bessel_y(T(3), eve::nan(eve::as<T>())), eve::nan(eve::as<T>()), 0);
  }

  for(int i=1; i < 5; ++i)
  {
    TTS_ULP_EQUAL(cyl_bessel_y(T(i), T(10)), T(std::cyl_neumann(i, 10)), 200);
    TTS_ULP_EQUAL(cyl_bessel_y(T(i), T(5)), T(std::cyl_neumann(i, 5)), 20);
    TTS_ULP_EQUAL(cyl_bessel_y(T(i), T(2)), T(std::cyl_neumann(i, 2)), 200);
    TTS_ULP_EQUAL(cyl_bessel_y(T(i), T(1)), T(std::cyl_neumann(i, 1)), 20);
    TTS_ULP_EQUAL(cyl_bessel_y(T(i), T(0)), T(std::cyl_neumann(i, 0)), 20);
    TTS_ULP_EQUAL(cyl_bessel_y(T(i), T(0.5)), T(std::cyl_neumann(i, 0.5)), 20);
  }

  if constexpr( eve::platform::supports_invalids )
  {
    TTS_IEEE_EQUAL(cyl_bessel_y(T(1), eve::nan(eve::as<T>()))  , eve::nan(eve::as<T>()) );
    TTS_IEEE_EQUAL(cyl_bessel_y(T(1), eve::inf(eve::as<T>()))  , T(0) );
  }

};
