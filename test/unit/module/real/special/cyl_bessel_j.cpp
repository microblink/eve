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
//#include <eve/function/diff/cyl_bessel_j.hpp>
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
// EVE_TEST_TYPES( "Check return types of cyl_bessel_j"
//             , eve::test::simd::ieee_reals
//             )
// <typename T>(eve::as<T>)
// {
//   using v_t = eve::element_type_t<T>;
//   using i_t = eve::as_integer_t<v_t>;
//   using I_t = eve::as_integer_t<T>;

//   TTS_EXPR_IS( eve::cyl_bessel_j(T(), T())  ,   T);
//   TTS_EXPR_IS( eve::cyl_bessel_j(v_t(),v_t()), v_t);
//   TTS_EXPR_IS( eve::cyl_bessel_j(i_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_j(I_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_j(I_t(),v_t()), T);
//   TTS_EXPR_IS( eve::cyl_bessel_j(i_t(),v_t()), v_t);
// };

//==================================================================================================
// cyl_bessel_j  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_j on wide"
        , eve::test::simd::ieee_doubles
        , eve::test::generate(eve::test::randoms(0.0, 10.0))
        )
<typename T>(T const&  )
{
  using v_t = eve::element_type_t<T>;
  auto eve__cyl_bessel_j =  [](v_t n, v_t x) {
    auto j = eve::cyl_bessel_j(n, x) ;
    return kumi::get<0>(j);
  };
  auto boost_cyl_bessel_j =  [](v_t n, v_t x) { return std::cyl_bessel_j(n, x); };
  if constexpr( eve::platform::supports_invalids )
  {
//     TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(2), eve::minf(eve::as<v_t>())), eve::zero(eve::as<v_t>()), 0);
//     TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(2), eve::inf(eve::as<v_t>())), eve::zero(eve::as<v_t>()), 0);
    TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(3), eve::nan(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
  }
  if constexpr(eve::cardinal_v<T> == 1)
  {
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(0), v_t(0)), v_t(boost_cyl_bessel_j(0.0, 0.0)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(1), v_t(0)), v_t(boost_cyl_bessel_j(1.0, 0.0)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(2), v_t(0)), v_t(boost_cyl_bessel_j(2.0, 0.0)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(0), v_t(1)), v_t(boost_cyl_bessel_j(0.0, 1.0)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(1), v_t(1)), v_t(boost_cyl_bessel_j(1.0, 1.0)), 10);
    for(int i=1; i < 2; i*= 2)
    {
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(i), v_t(10)), v_t(boost_cyl_bessel_j(i, 10)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(i), v_t(5)), v_t(boost_cyl_bessel_j(i, 5)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(i), v_t(2)), v_t(boost_cyl_bessel_j(i, 2)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(i), v_t(1)), v_t(boost_cyl_bessel_j(i, 1)), 10);
      TTS_ULP_EQUAL(eve__cyl_bessel_j(v_t(i), v_t(0)), v_t(boost_cyl_bessel_j(i, 0)), 10);
    }
  }
//  using v_t = eve::element_type_t<T>;
//   using eve::cyl_bessel_j;
//   using eve::as;
//   for(int i=1; i < 4 ; ++i)
//   {
//     TTS_ULP_EQUAL( cyl_bessel_j(i, a0),  map([i](auto e){return boost::math::cyl_bessel_j(i, e);}, a0), 5);
//     auto dcyl_bessel_j = [i](auto e){return v_t( -boost::math::cyl_bessel_j(i-1, e));};
//     TTS_ULP_EQUAL( eve::diff(cyl_bessel_j)(i, a0),  map(dcyl_bessel_j, a0), 5);
//   }

//   if constexpr( eve::platform::supports_invalids )
//   {
//     TTS_IEEE_EQUAL(cyl_bessel_j(T(1), eve::nan(eve::as<T>()))  , eve::nan(eve::as<T>()) );
//     TTS_IEEE_EQUAL(cyl_bessel_j(T(1), eve::inf(eve::as<T>()))   , T(0) );
//   }


//   for(int i=1; i < 4 ; ++i)
//   {
//     TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(0))  , eve::rec(T(i-1)), 0.5);
//     TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(0.5)), T(boost::math::cyl_bessel_j(i, 0.5)), 2.0);
//     TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(1))  , T(boost::math::cyl_bessel_j(i, 1.0)), 4.0);
//     TTS_ULP_EQUAL(cyl_bessel_j(T(i), T(10)) , T(boost::math::cyl_bessel_j(i, 10.0)), 0.5);
//   }
//   for(int i=1; i < 4 ; ++i)
//   {
//     TTS_ULP_EQUAL(cyl_bessel_j(i, T(0))  , eve::rec(T(i-1)), 0.5);
//     TTS_ULP_EQUAL(cyl_bessel_j(i, T(0.5)), T(boost::math::cyl_bessel_j(i, 0.5)), 2.0);
//     TTS_ULP_EQUAL(cyl_bessel_j(i, T(1))  , T(boost::math::cyl_bessel_j(i, 1.0)), 4.0);
//     TTS_ULP_EQUAL(cyl_bessel_j(i, T(10)) , T(boost::math::cyl_bessel_j(i, 10.0)), 0.5);
//   }
//   using elt_t =  eve::element_type_t<T>;

//   TTS_ULP_EQUAL(cyl_bessel_j(elt_t(2.0), elt_t(0.5)), (boost::math::cyl_bessel_j(elt_t(2), elt_t(0.5))), 2.0);
//   TTS_ULP_EQUAL(cyl_bessel_j(elt_t(6000), elt_t(0.5)), (boost::math::cyl_bessel_j(elt_t(6000), elt_t(0.5))), 3.0);
};
