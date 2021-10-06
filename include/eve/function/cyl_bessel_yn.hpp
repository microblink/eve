//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/detail/overload.hpp>

namespace eve
{
  //================================================================================================
  //! @addtogroup special
  //! @{
  //! @var cyl_bessel_yn
  //!
  //! @brief Callable object computing the bessel regular function of the second kind:\f$\mbox{J}_\nu(x)\f$, \f$\mbox{Y}_\nu(x)\f$,
  //! and their first order derivatives for integral values of nu.
  //!
  //! **Required header:** `#include <eve/function/cyl_bessel_yn.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the cyl_bessel_yn operation                                |
  //!
  //! ---
  //!
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
  //!  auto operator()( real_value auto nu, floating_real_value auto x ) const noexcept
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //!
  //! **Parameters**
  //!
  //!`nu`:   [real value](@ref eve::real_value). nu must be a [flint}(@ref eve::is_flint), but can be negative
  //!
  //!`x`:   [floating real-value](@ref eve::floating_real_value). `x` must be positive.
  //!
  //! **Return value**
  //!
  //! Returns [elementwise](@ref glossary_elementwise) thhe value of the Neumann bessel function of integral oreder nu.
  //!
  //!  \f$\displaystyle \mbox{Y}_\nu(x) = \frac{\mbox{J}_\nu(x)\cos(\pi x)-\mbox{J}_{-\nu(x)}}{sin(\nu\pi}\f$.
  //!
  //! These functions return 0 for infinite positive `x`
  //!
  //! For big values of nu` or/and `x` the computation may fail to converge. In this case nan is returned.
  //!
  //! ---
  //!
  //! #### Example
  //!
  //! @godbolt{doc/core/cyl_bessel_yn.cpp}
  //!
  //!  @}
  //================================================================================================
  EVE_MAKE_CALLABLE(cyl_bessel_yn_, cyl_bessel_yn);
}

#include <eve/module/real/special/function/regular/generic/cyl_bessel_yn.hpp>
