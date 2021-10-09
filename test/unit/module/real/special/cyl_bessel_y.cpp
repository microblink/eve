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
// EVE_TEST_TYPES( "Check return types of cyl_bessel_y"
//             , eve::test::simd::ieee_reals
//             )
// <typename T>(eve::as<T>)
// {
//   using v_t = eve::element_type_t<T>;
//   using i_t = eve::as_integer_t<v_t>;
//   using I_t = eve::as_integer_t<T>;

//   TTS_EXPR_IS( eve::cyl_bessel_y(T(), T())  ,   T);
//   TTS_EXPR_IS( eve::cyl_bessel_y(v_t(),v_t()), v_t);
//   TTS_EXPR_IS( eve::cyl_bessel_y(i_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_y(I_t(),T()),   T);
//   TTS_EXPR_IS( eve::cyl_bessel_y(I_t(),v_t()), T);
//   TTS_EXPR_IS( eve::cyl_bessel_y(i_t(),v_t()), v_t);
// };

//==================================================================================================
//=== cyl_bessel_y  tests
//==================================================================================================
EVE_TEST( "Check behavior of cyl_bessel_y on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(0.0, 10.0),
                              eve::test::ramp(0))
        )
<typename T>(T a0, T )
{
  using v_t = eve::element_type_t<T>;
  auto eve__cyl_bessel_y =  [](auto n, auto x) {
    return eve::cyl_bessel_y(n, x) ;
  };
  auto boost_cyl_bessel_y =  [](v_t n, v_t x) { return std::cyl_neumann(n, x); };
//   if constexpr( eve::platform::supports_invalids )
//   {
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(2), eve::minf(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(2), eve::inf(eve::as<v_t>())), eve::zero(eve::as<v_t>()), 0);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(3), eve::nan(eve::as<v_t>())), eve::nan(eve::as<v_t>()), 0);

//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), eve::minf(eve::as<T>())), eve::nan(eve::as<T>()), 0);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), eve::inf(eve::as<T>())), eve::zero(eve::as<T>()), 0);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(3), eve::nan(eve::as<T>())), eve::nan(eve::as<T>()), 0);
//   }
//   if constexpr(eve::cardinal_v<T> == 1)
//   {
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(0), v_t(0)), v_t(boost_cyl_bessel_y(v_t(0.0), v_t(0.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(1), v_t(0)), v_t(boost_cyl_bessel_y(v_t(1.0), v_t(0.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(2), v_t(0)), v_t(boost_cyl_bessel_y(v_t(2.0), v_t(0.0))), 10);

//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), T(0)), T(boost_cyl_bessel_y(v_t(0.0), v_t(0.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), T(0)), T(boost_cyl_bessel_y(v_t(1.0), v_t(0.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), T(0)), T(boost_cyl_bessel_y(v_t(2.0), v_t(0.0))), 10);

//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(0.1), v_t(1)), v_t(boost_cyl_bessel_y(v_t(0.1), v_t(1.0))), 300);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(0), v_t(10)), v_t(boost_cyl_bessel_y(v_t(0.0), v_t(10.0))), 30);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(0), v_t(1)), v_t(boost_cyl_bessel_y(v_t(0.0), v_t(1.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(1), v_t(1)), v_t(boost_cyl_bessel_y(v_t(1.0), v_t(1.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(2), v_t(1)), v_t(boost_cyl_bessel_y(v_t(2.0), v_t(1.0))), 10);

//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), T(1)), T(boost_cyl_bessel_y(v_t(0.0), v_t(1.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), T(1)), T(boost_cyl_bessel_y(v_t(1.0), v_t(1.0))), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), T(1)), T(boost_cyl_bessel_y(v_t(2.0), v_t(1.0))), 10);

//     for(int ii=1; ii < 8; ii*= 2)
//     {
//       auto i = ii+v_t(0.5);
//       std::cout << "i " << i << std::endl;
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(10)), v_t(boost_cyl_bessel_y(i, 10)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(5)), v_t(boost_cyl_bessel_y(i, 5)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(2)), v_t(boost_cyl_bessel_y(i, 2)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(1)), v_t(boost_cyl_bessel_y(i, 1)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(0)), v_t(boost_cyl_bessel_y(i, 0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), eve::eps(eve::as<v_t>())), v_t(boost_cyl_bessel_y(i, eve::eps(eve::as<v_t>()))), 10);

//       TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(10)), T(boost_cyl_bessel_y(i, 10)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(5)), T(boost_cyl_bessel_y(i, 5)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(2)), T(boost_cyl_bessel_y(i, 2)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(1)), T(boost_cyl_bessel_y(i, 1)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(0)), T(boost_cyl_bessel_y(i, 0)), 10);
//       TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), eve::eps(eve::as<T>())), T(boost_cyl_bessel_y(i, eve::eps(eve::as<v_t>()))), 10);
//     }
//   }
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), T(0)), T(boost_cyl_bessel_y(0.0, 0.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), T(0)), T(boost_cyl_bessel_y(1.0, 0.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), T(0)), T(boost_cyl_bessel_y(2.0, 0.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), T(1)), T(boost_cyl_bessel_y(0.0, 1.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), T(1)), T(boost_cyl_bessel_y(1.0, 1.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), T(1)), T(boost_cyl_bessel_y(2.0, 1.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), T(3)), T(boost_cyl_bessel_y(0.0, 3.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), T(3)), T(boost_cyl_bessel_y(1.0, 3.0)), 10);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), T(3)), T(boost_cyl_bessel_y(2.0, 3.0)), 10);

//  a0/= 100;
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), a0), map(boost_cyl_bessel_y, T(0), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), a0), map(boost_cyl_bessel_y, T(1), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), a0), map(boost_cyl_bessel_y, T(2), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0.5), a0), map(boost_cyl_bessel_y, T(0.5), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(n, a0), map(boost_cyl_bessel_y, n, a0), 30);
  a0+= 500;
//  TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), a0), map(boost_cyl_bessel_y, T(0), a0), 30);
  TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), a0), map(boost_cyl_bessel_y, T(1), a0), 30);
//  TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), a0), map(boost_cyl_bessel_y, T(2), a0), 30);
//  TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0.5), a0), map(boost_cyl_bessel_y, T(0.5), a0), 30);
//  TTS_ULP_EQUAL(eve__cyl_bessel_y(n, a0), map(boost_cyl_bessel_y, n, a0), 30);
//   a0+= 15000;
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0), a0), map(boost_cyl_bessel_y, T(0), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(1), a0), map(boost_cyl_bessel_y, T(1), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(2), a0), map(boost_cyl_bessel_y, T(2), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(T(0.5), a0), map(boost_cyl_bessel_y, T(0.5), a0), 30);
//   TTS_ULP_EQUAL(eve__cyl_bessel_y(n, a0), map(boost_cyl_bessel_y, n, a0), 30);
};

// EVE_TEST_TYPES( "Check return types of cyl_bessel_y"
//             , eve::test::simd::ieee_doubles
//             )
// <typename T>(eve::as<T>)
// {
//   using v_t = eve::element_type_t<T>;
//   auto eve__cyl_bessel_y =  [](auto n, auto x) { return eve::cyl_bessel_y(n, x);  };
//   auto boost_cyl_bessel_y = [](v_t n, v_t x) -> v_t{ return std::cyl_neumann(n, x); };

// //  TTS_ULP_EQUAL(eve__cyl_bessel_y(1.5, 2.0), v_t(boost_cyl_bessel_y(1.5, 2.0)), 10);
//   for(int ii=1; ii < 8; ii*= 2)
//   {
//     auto i = ii+v_t(0.5);
//     std::cout << "i " << i << std::endl;
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(10)), v_t(boost_cyl_bessel_y(i, 10)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(5)), v_t(boost_cyl_bessel_y(i, 5)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(2)), v_t(boost_cyl_bessel_y(i, 2)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(1)), v_t(boost_cyl_bessel_y(i, 1)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), v_t(0)), v_t(boost_cyl_bessel_y(i, 0)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(v_t(i), eve::eps(eve::as<v_t>())), v_t(boost_cyl_bessel_y(i, eve::eps(eve::as<v_t>()))), 10);

//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(10)), T(boost_cyl_bessel_y(i, 10)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(5)), T(boost_cyl_bessel_y(i, 5)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(2)), T(boost_cyl_bessel_y(i, 2)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(1)), T(boost_cyl_bessel_y(i, 1)), 10);
//     TTS_ULP_EQUAL(eve__cyl_bessel_y(T(i), T(0)), T(boost_cyl_bessel_y(i, 0)), 10);
//   }
// };
