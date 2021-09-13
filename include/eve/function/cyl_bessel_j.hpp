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
  //! @var cyl_bessel_j
  //!
  //! @brief Callable object computing the bessel regular funcyion of the first  kind:\f$\mbox{J}_\nu(x)\f$,.
  //!
  //! **Required header:** `#include <eve/function/cyl_bessel_j.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the cyl_bessel_j operation                                 |
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
  //!`x`:   [floating real-value](@ref eve::floating_real_value). `x` must be positive if `nu` is not a [flint](@ref eve::is_flint).
  //!
  //! **Return value**
  //!
  //! Returns [elementwise](@ref glossary_elementwise)  \f$\mbox{J}_\nu(x)\f$.
  //!
  //! With \f$\displaystyle \mbox{J}_\nu(x) = \left(\frac{z}2\right)^\nu \sum_{k = 0}^\infty \frac{(-z^2/4)^k}{k!\Gamma(\nu+k+1)}\f$
  //!
  //! ---
  //!
  //! #### Supported decorators
  //!
 //!   * eve::diff, eve::diff_1st, eve::diff_nth
  //!
  //!     **Required header:** `#include <eve/function/diff/cyl_bessel_j.hpp>`
  //!
  //!     The expression `eve::diff(eve::cyl_bessel_j)(x)` computes the derivative of the function at `x`.
  //!
  //! ---
  //!
  //! #### Example
  //!
  //! @godbolt{doc/core/cyl_bessel_j.cpp}
  //!
  //!  @}
  //================================================================================================
   EVE_MAKE_CALLABLE(cyl_bessel_j_, cyl_bessel_j);
}

#include <eve/module/real/special/function/regular/generic/cyl_bessel_j.hpp>
