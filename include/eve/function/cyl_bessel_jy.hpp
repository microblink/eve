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
  //! @var cyl_bessel_jy
  //!
  //! @brief Callable object computing the bessel regular funcyion of the first and second kind:\f$\mbox{J}_\nu(x)\f$, \f$\mbox{Y}_\nu(x)\f$,
  //! and their first order derivatives.
  //!
  //! **Required header:** `#include <eve/function/cyl_bessel_jy.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the cyl_bessel_jy operation                                |
  //!
  //! ---
  //!
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
  //!  auto operator()( real_value auto nu, floating_real_value auto x ) const noexcept
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //!
  //! **Parameters**
  //!
  //!`nu`:   [real value](@ref eve::real_value).
  //!
  //!`x`:   [floating real-value](@ref eve::floating_real_value). `x` must be positive.
  //!
  //!@warning
  //! On the real field f$\mbox{J}_\nu(x)\f$ is the only bessel function defined for negative inputs and only if nu is a [flint](@ref eve::is_flint).
  //! In this case, a direct call to `cyl_bessel_j(nu, x)` will return the expected value i.e. \f$(-1)^nu  \mbox{J}_\nu(-x)\f$, all the  outputs of
  //! `cyl_bessel_jy(nu, x)` will be nan.
  //!
  //! **Return value**
  //!
  //! Returns [elementwise](@ref glossary_elementwise) a kumi quadruplet consisting of the respective values:
  //! \f$\mbox{J}_\nu(x)\f$, \f$\mbox{J}_\nu'(x)\f$, \f$\mbox{Y}_\nu(x)\f$, \f$\mbox{Y}_\nu'(x)\f$.
  //!
  //! With \f$\displaystyle \mbox{J}_\nu(x) = \left(\frac{z}2\right)^\nu \sum_{k = 0}^\infty \frac{(-z^2/4)^k}{k!\Gamma(\nu+k+1)}\f$
  //!  and \f$\displaystyle \mbox{Y}_\nu(x) = \frac{\mbox{J}_\nu(x)\cos(\pi x)-\mbox{J}_{-\nu(x)}}{sin(\nu\pi}\f$.
  //!
  //! These functions return 0 for infinite positive `x`
  //!
  //! For big values of nu` or/and `x` the computation may fail to converge. In this case nan is returned.
  //!
  //! ---
  //!
  //! #### Example
  //!
  //! @godbolt{doc/core/cyl_bessel_jy.cpp}
  //!
  //!  @}
  //================================================================================================
  EVE_MAKE_CALLABLE(cyl_bessel_jy_, cyl_bessel_jy);
}

#include <eve/module/real/special/function/regular/generic/cyl_bessel_jy.hpp>
