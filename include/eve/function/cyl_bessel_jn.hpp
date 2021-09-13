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
  //! @var cyl_bessel_jn
  //!
  //! @brief Callable object computing the bessel regular funcyion of the first and second kind:\f$\mbox{J}_\nu(x)\f$, \f$\mbox{N}_\nu(x)\f$,
  //! and their first order derivatives.
  //!
  //! **Required header:** `#include <eve/function/cyl_bessel_jn.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the cyl_bessel_jn operation                                |
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
  //! On the real field f$\mbox{J}_\nu(x)\f$ is the only bessel function defined for negative inputs and so only if nu is a [flint](@ref eve::is_flint).
  //! In this case, a direct call to `cyl_bessel_j(nu, x)` will return the expected value i.e. \f$(-1)^nu  \mbox{J}_\nu(-x)\f$, but the first output of
  //! `cyl_bessel_jn(nu, x)` will be a nan.
  //!
  //! **Return value**
  //!
  //! Returns [elementwise](@ref glossary_elementwise) a kumi quadruplet consisting of the respective values:
  //! \f$\mbox{J}_\nu(x)\f$, \f$\mbox{J}_\nu'(x)\f$, \f$\mbox{N}_\nu(x)\f$, \f$\mbox{N}_\nu'(x)\f$.
  //!
  //! With \f$\displaystyle \mbox{J}_\nu(x) = \left(\frac{z}2\right)^\nu \sum_{k = 0}^\infty \frac{(-z^2/4)^k}{k!\Gamma(\nu+k+1)}\f$
  //!  and \f$\displaystyle \mbox{N}_\nu(x) = \frac{\mbox{J}_\nu(x)\cos(\pi x)-\mbox{J}_{-\nu(x)}}{sin(\nu\pi}\f$.
  //!
  //! ---
  //!
  //! #### Example
  //!
  //! @godbolt{doc/core/cyl_bessel_jn.cpp}
  //!
  //!  @}
  //================================================================================================
  EVE_MAKE_CALLABLE(cyl_bessel_jn_, cyl_bessel_jn);
}

#include <eve/module/real/special/function/regular/generic/cyl_bessel_jn.hpp>
