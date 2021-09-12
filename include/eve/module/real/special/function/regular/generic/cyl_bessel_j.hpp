//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/concept/value.hpp>
#include <eve/function/fms.hpp>
#include <eve/function/abs.hpp>
#include <eve/function/any.hpp>
#include <eve/function/cospi.hpp>
#include <eve/function/cyl_bessel_j0.hpp>
#include <eve/function/cyl_bessel_j1.hpp>
#include <eve/detail/hz_device.hpp>
#include <eve/function/fms.hpp>
#include <eve/function/fnma.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/is_ltz.hpp>
#include <eve/function/is_nltz.hpp>
#include <eve/function/is_nlez.hpp>
#include <eve/function/is_flint.hpp>
#include <eve/function/is_odd.hpp>
#include <eve/constant/true.hpp>
#include <type_traits>

#include <eve/constant/eps.hpp>
#include <eve/constant/inf.hpp>
#include <eve/constant/invpi.hpp>
#include <eve/constant/nan.hpp>
#include <eve/constant/pi.hpp>
#include <eve/constant/smallestposval.hpp>
#include <eve/function/abs.hpp>
#include <eve/function/average.hpp>
#include <eve/function/copysign.hpp>
#include <eve/function/fam.hpp>
#include <eve/function/fma.hpp>
#include <eve/function/is_eqz.hpp>
#include <eve/function/log.hpp>
#include <eve/function/rec.hpp>
#include <eve/function/sinpic.hpp>
#include <eve/function/sinhc.hpp>
#include <eve/function/sqr.hpp>
#include <eve/function/sqrt.hpp>
#include <eve/function/tgamma.hpp>

#include <eve/function/sin.hpp>
#include <eve/function/cosh.hpp>

namespace eve::detail
{

  template<real_scalar_value I, floating_real_scalar_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I n, T x) noexcept
  {
    if constexpr(std::is_integral_v<I>) return cyl_bessel_j(T(n), x);
    else
    {
      T nu(n);
      T jnu, jpnu, nnu, npnu;

      if (is_eqz(x))
      {
        if (is_eqz(nu))
        {
          return kumi::make_tuple(T(1), T(0), inf(as(x)), inf(as(x))); //jnu, jpnu, nnu, npnu
        }
        else if (nu ==one(as(x)))
        {
          return kumi::make_tuple(T(0), T(0.5), inf(as(x)), inf(as(x))); //jnu, jpnu, nnu, npnu
        }
        else
        {
          return kumi::make_tuple(T(0), T(0), inf(as(x)), inf(as(x))); //jnu, jpnu, nnu, npnu
        }
      }
      const T eps = eve::eps(as(x));
      //  When the multiplier is n i.e.
      //  fp_min = n * min()
      //  Then j_0 and n_0 tank at x = 8 * n (j_0 = 0 and n_0 = nan)!
      //const T fp_min = T(20) * std::numeric_limits<T>::min();

      const T fp_min = eve::sqrt(smallestposval(as(x)));
      constexpr int max_iter = 15000;
      const T x_min = T(2);
      const int nl = (x < x_min
                      ? static_cast<int>(nu + T(0.5))
                      : eve::max(0, static_cast<int>(nu - x + T(1.5))));
      const T mu = nu - nl;
      const T mu2 = sqr(mu);
      const T xi = rec(x);
      const T xi2 = xi + xi;
      T w = xi2 *invpi(as(x));
      int isign = 1;
      T h = nu * xi;
      if (h < fp_min)  h = fp_min;
      T b = xi2 * nu;
      T d = T(0);
      T c = h;
      int i;
      for (i = 1; i <= max_iter; ++i)
        {
          b += xi2;
          d = b - d;
          if (std::abs(d) < fp_min) d = fp_min;
          c = b - T(1) / c;
          if (std::abs(c) < fp_min) c = fp_min;
          d = T(1) / d;
          const T del = c * d;
          h *= del;
          if (d < T(0))  isign = -isign;
          if (eve::abs(del - T(1)) < eve::eps(as(x))) break;
        }
//      std::cout << "i " << i << std::endl;
//     if (i > max_iter) return nan(as(x));
//         std::throw_runtime_error(n("Argument x too large in bessel_jn; "
//                                        "try asymptotic expansion."));
//       std::cout << "isign "<< isign << std::endl;
//       std::cout << "fp_min  "<< fp_min  << std::endl;
      T jnul = isign * fp_min;
      T jpnul = h * jnul;
      T jnul1 = jnul;
      T jpnu1 = jpnul;
      T fact = nu * xi;
//       std::cout << "a jpnul "<< jpnul << std::endl;
//       std::cout << "a jnul  "<< jnul  << std::endl;
      for ( int l = nl; l >= 1; --l )
        {
          const T jnutemp = fma(fact, jnul, jpnul);
          fact -= xi;
          jpnul = fact * jnutemp - jnul;
          jnul = jnutemp;
        }
      if (jnul == T(0)) jnul = eve::eps(as(x));
      T f = jpnul / jnul;
//       std::cout << "jpnul "<< jpnul << std::endl;
//       std::cout << "jnul  "<< jnul  << std::endl;
      T nmu, nnu1, npmu, jmu;
      if (x < x_min)
      {
//        std::cout << "icitte" << std::endl;
         const T x2 = x / T(2);
         const T pimu = pi(as(x)) * mu;
          T fact = (std::abs(pimu) < eps
                      ? T(1) : pimu / std::sin(pimu));
          T d = -std::log(x2);
          T e = mu * d;
          T fact2 = (std::abs(e) < eps
                       ? T(1) : std::sinh(e) / e);
          T gam1, gam2, gampl, gammi;

          auto gamma_temme = [&gam1, &gam2, &gampl, &gammi](auto mu) {
            auto gamma_e = T(0.57721566490153286060651209008240243104215933593992);
            gampl = rec(tgamma(inc(mu)));
            gammi = rec(tgamma(oneminus(mu)));
            gam1 = (eve::abs(mu) < eve::eps(as(mu))) ? gamma_e : (gammi - gampl) /(mu+mu);
            gam2 = average(gammi, gampl);
            return;
          };

          gamma_temme(mu);
//          std::cout << mu << " -- " <<  gam1<< " -- " << gam2<< " -- " << gampl<< " -- " << gammi << std::endl;
          //gamma_temme(mu, gam1, gam2, gampl, gammi);
          T ff = (T(2) / pi(as(x)))
            * fact * (gam1 * eve::cosh(e) + gam2 * fact2 * d);
          e = std::exp(e);
          T p = e / (pi(as(x)) * gampl);
          T q = T(1) / (e * pi(as(x)) * gammi);
          const T pimu2 = pimu / T(2);
          T fact3 = (std::abs(pimu2) < eps
                     ? T(1) : eve::sin(pimu2) / pimu2 );
          T r = pi(as(x)) * pimu2 * fact3 * fact3;
          T c = T(1);
          d = -x2 * x2;
          T sum = ff + r * q;
          T sum1 = p;
          for (i = 1; i <= max_iter; ++i)
          {
            ff = (i * ff + p + q) / (i * i - mu2);
            c *= d / T(i);
            p /= T(i) - mu;
            q /= T(i) + mu;
            const T del = c * (ff + r * q);
            sum += del;
            const T del1 = c * p - i * del;
            sum1 += del1;
            if ( eve::abs(del) < eve::eps(as(x)) * (T(1) + eve::abs(sum)) )
              break;
          }
//          std::cout << "2 i "<< i << std::endl;
//           if ( i > max_iter )
//             std::throw_runtime_error(n("Bessel y series failed to converge "
//                                            "in bessel_jn."));
          nmu = -sum;
          nnu1 = -sum1 * xi2;
          npmu = mu * xi * nmu - nnu1;
          jmu = w / (npmu - f * nmu);
//         std::cout << "sum   "<< sum << std::endl;
//         std::cout << "sum1  "<< sum1<< std::endl;
//         std::cout << "nnu1  "<< nnu1<< std::endl;
//         std::cout << "nmu   "<< nmu << std::endl;
//         std::cout << "jmu   "<< jmu << std::endl;
//         std::cout << "jw    "<< w   << std::endl;
//         std::cout << "npmu  "<< npmu<< std::endl;
//         std::cout << "f     "<< f   << std::endl;


      }
      else
      {
//        std::cout << "latte" << std::endl;
        T a = T(0.25L) - mu2;
        T q = T(1);
        T p = -xi / T(2);
        T br = T(2) * x;
        T bi = T(2);
        T fact = a * xi / (p * p + q * q);
        T cr = br + q * fact;
        T ci = bi + p * fact;
        T den = br * br + bi * bi;
        T dr = br / den;
        T di = -bi / den;
        T dlr = cr * dr - ci * di;
        T dli = cr * di + ci * dr;
        T temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        int i;
        for (i = 2; i <= max_iter; ++i)
        {
          a += T(2 * (i - 1));
          bi += T(2);
          dr = a * dr + br;
          di = a * di + bi;
          if (eve::abs(dr) + eve::abs(di) < fp_min) dr = fp_min;
          fact = a / (cr * cr + ci * ci);
          cr = br + cr * fact;
          ci = bi - ci * fact;
          if (eve::abs(cr) + eve::abs(ci) < fp_min) cr = fp_min;
          den = dr * dr + di * di;
          dr /= den;
          di /= -den;
          dlr = cr * dr - ci * di;
          dli = cr * di + ci * dr;
          temp = p * dlr - q * dli;
          q = p * dli + q * dlr;
          p = temp;
          if (eve::abs(dlr - T(1)) + eve::abs(dli) < eps)  break;
        }
//        if (i > max_iter) return nan(as(x));
//             std::throw_runtime_error(n("Lentz's method failed "
//                                            "in bessel_jn."));
        const T gam = (p - f) / q;
        jmu = std::sqrt(w / ((p - f) * gam + q));
//        std::cout << "jmu   "<< jmu << std::endl;
        jmu = eve::copysign(jmu, jnul);
        nmu = gam * jmu;
        npmu = (p + q / gam) * nmu;
        nnu1 = mu * xi * nmu - npmu;
      }
      fact = jmu / jnul;
      jnu = fact * jnul1;
 //          std::cout << "jnu   "<< jnu << std::endl;
//           std::cout << "jnul1 "<< jnul1<< std::endl;
//           std::cout << "fact "<< fact << std::endl;

      jpnu = fact * jpnu1;
      for (i = 1; i <= nl; ++i)
      {
        const T nnutemp = (mu + i) * xi2 * nnu1 - nmu;
        nmu = nnu1;
        nnu1 = nnutemp;
      }
      nnu = nmu;
      npnu = nu * xi * nmu - nnu1;
      return kumi::make_tuple(jnu, jpnu, nnu, npnu);
    }
  }


  template<floating_real_simd_value T>
  EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), T nu, T x) noexcept
  {
    T jnu, jpnu, nnu, npnu;
    
      if (is_eqz(x))
      {
        if (is_eqz(nu))
        {
          return kumi::make_tuple(T(1), T(0), inf(as(x)), inf(as(x))); //jnu, jpnu, nnu, npnu
        }
        else if (nu ==one(as(x)))
        {
          return kumi::make_tuple(T(0), T(0.5), inf(as(x)), inf(as(x))); //jnu, jpnu, nnu, npnu
        }
        else
        {
          return kumi::make_tuple(T(0), T(0), inf(as(x)), inf(as(x))); //jnu, jpnu, nnu, npnu
        }
      }
      const T eps = eve::eps(as(x));
      //  When the multiplier is n i.e.
      //  fp_min = n * min()
      //  Then j_0 and n_0 tank at x = 8 * n (j_0 = 0 and n_0 = nan)!
      //const T fp_min = T(20) * std::numeric_limits<T>::min();

      const T fp_min = eve::sqrt(smallestposval(as(x)));
      constexpr int max_iter = 15000;
      const T x_min = T(2);
      const int nl = (x < x_min
                      ? static_cast<int>(nu + T(0.5))
                      : eve::max(0, static_cast<int>(nu - x + T(1.5))));
      const T mu = nu - nl;
      const T mu2 = sqr(mu);
      const T xi = rec(x);
      const T xi2 = xi + xi;
      T w = xi2 *invpi(as(x));
      int isign = 1;
      T h = nu * xi;
      if (h < fp_min)  h = fp_min;
      T b = xi2 * nu;
      T d = T(0);
      T c = h;
      int i;
      for (i = 1; i <= max_iter; ++i)
        {
          b += xi2;
          d = b - d;
          if (std::abs(d) < fp_min) d = fp_min;
          c = b - T(1) / c;
          if (std::abs(c) < fp_min) c = fp_min;
          d = T(1) / d;
          const T del = c * d;
          h *= del;
          if (d < T(0))  isign = -isign;
          if (eve::abs(del - T(1)) < eve::eps(as(x))) break;
        }
//      std::cout << "i " << i << std::endl;
//     if (i > max_iter) return nan(as(x));
//         std::throw_runtime_error(n("Argument x too large in bessel_jn; "
//                                        "try asymptotic expansion."));
//       std::cout << "isign "<< isign << std::endl;
//       std::cout << "fp_min  "<< fp_min  << std::endl;
      T jnul = isign * fp_min;
      T jpnul = h * jnul;
      T jnul1 = jnul;
      T jpnu1 = jpnul;
      T fact = nu * xi;
//       std::cout << "a jpnul "<< jpnul << std::endl;
//       std::cout << "a jnul  "<< jnul  << std::endl;
      for ( int l = nl; l >= 1; --l )
        {
          const T jnutemp = fma(fact, jnul, jpnul);
          fact -= xi;
          jpnul = fact * jnutemp - jnul;
          jnul = jnutemp;
        }
      if (jnul == T(0)) jnul = eve::eps(as(x));
      T f = jpnul / jnul;
//       std::cout << "jpnul "<< jpnul << std::endl;
//       std::cout << "jnul  "<< jnul  << std::endl;
      T nmu, nnu1, npmu, jmu;
      if (x < x_min)
      {
//        std::cout << "icitte" << std::endl;
         const T x2 = x / T(2);
         const T pimu = pi(as(x)) * mu;
          T fact = (std::abs(pimu) < eps
                      ? T(1) : pimu / std::sin(pimu));
          T d = -std::log(x2);
          T e = mu * d;
          T fact2 = (std::abs(e) < eps
                       ? T(1) : std::sinh(e) / e);
          T gam1, gam2, gampl, gammi;

          auto gamma_temme = [&gam1, &gam2, &gampl, &gammi](auto mu) {
            auto gamma_e = T(0.57721566490153286060651209008240243104215933593992);
            gampl = rec(tgamma(inc(mu)));
            gammi = rec(tgamma(oneminus(mu)));
            gam1 = (eve::abs(mu) < eve::eps(as(mu))) ? gamma_e : (gammi - gampl) /(mu+mu);
            gam2 = average(gammi, gampl);
            return;
          };

          gamma_temme(mu);
//          std::cout << mu << " -- " <<  gam1<< " -- " << gam2<< " -- " << gampl<< " -- " << gammi << std::endl;
          //gamma_temme(mu, gam1, gam2, gampl, gammi);
          T ff = (T(2) / pi(as(x)))
            * fact * (gam1 * eve::cosh(e) + gam2 * fact2 * d);
          e = std::exp(e);
          T p = e / (pi(as(x)) * gampl);
          T q = T(1) / (e * pi(as(x)) * gammi);
          const T pimu2 = pimu / T(2);
          T fact3 = (std::abs(pimu2) < eps
                     ? T(1) : eve::sin(pimu2) / pimu2 );
          T r = pi(as(x)) * pimu2 * fact3 * fact3;
          T c = T(1);
          d = -x2 * x2;
          T sum = ff + r * q;
          T sum1 = p;
          for (i = 1; i <= max_iter; ++i)
          {
            ff = (i * ff + p + q) / (i * i - mu2);
            c *= d / T(i);
            p /= T(i) - mu;
            q /= T(i) + mu;
            const T del = c * (ff + r * q);
            sum += del;
            const T del1 = c * p - i * del;
            sum1 += del1;
            if ( eve::abs(del) < eve::eps(as(x)) * (T(1) + eve::abs(sum)) )
              break;
          }
//          std::cout << "2 i "<< i << std::endl;
//           if ( i > max_iter )
//             std::throw_runtime_error(n("Bessel y series failed to converge "
//                                            "in bessel_jn."));
          nmu = -sum;
          nnu1 = -sum1 * xi2;
          npmu = mu * xi * nmu - nnu1;
          jmu = w / (npmu - f * nmu);
//         std::cout << "sum   "<< sum << std::endl;
//         std::cout << "sum1  "<< sum1<< std::endl;
//         std::cout << "nnu1  "<< nnu1<< std::endl;
//         std::cout << "nmu   "<< nmu << std::endl;
//         std::cout << "jmu   "<< jmu << std::endl;
//         std::cout << "jw    "<< w   << std::endl;
//         std::cout << "npmu  "<< npmu<< std::endl;
//         std::cout << "f     "<< f   << std::endl;


      }
      else
      {
//        std::cout << "latte" << std::endl;
        T a = T(0.25L) - mu2;
        T q = T(1);
        T p = -xi / T(2);
        T br = T(2) * x;
        T bi = T(2);
        T fact = a * xi / (p * p + q * q);
        T cr = br + q * fact;
        T ci = bi + p * fact;
        T den = br * br + bi * bi;
        T dr = br / den;
        T di = -bi / den;
        T dlr = cr * dr - ci * di;
        T dli = cr * di + ci * dr;
        T temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        int i;
        for (i = 2; i <= max_iter; ++i)
        {
          a += T(2 * (i - 1));
          bi += T(2);
          dr = a * dr + br;
          di = a * di + bi;
          if (eve::abs(dr) + eve::abs(di) < fp_min) dr = fp_min;
          fact = a / (cr * cr + ci * ci);
          cr = br + cr * fact;
          ci = bi - ci * fact;
          if (eve::abs(cr) + eve::abs(ci) < fp_min) cr = fp_min;
          den = dr * dr + di * di;
          dr /= den;
          di /= -den;
          dlr = cr * dr - ci * di;
          dli = cr * di + ci * dr;
          temp = p * dlr - q * dli;
          q = p * dli + q * dlr;
          p = temp;
          if (eve::abs(dlr - T(1)) + eve::abs(dli) < eps)  break;
        }
//        if (i > max_iter) return nan(as(x));
//             std::throw_runtime_error(n("Lentz's method failed "
//                                            "in bessel_jn."));
        const T gam = (p - f) / q;
        jmu = std::sqrt(w / ((p - f) * gam + q));
//        std::cout << "jmu   "<< jmu << std::endl;
        jmu = eve::copysign(jmu, jnul);
        nmu = gam * jmu;
        npmu = (p + q / gam) * nmu;
        nnu1 = mu * xi * nmu - npmu;
      }
      fact = jmu / jnul;
      jnu = fact * jnul1;
 //          std::cout << "jnu   "<< jnu << std::endl;
//           std::cout << "jnul1 "<< jnul1<< std::endl;
//           std::cout << "fact "<< fact << std::endl;

      jpnu = fact * jpnu1;
      for (i = 1; i <= nl; ++i)
      {
        const T nnutemp = (mu + i) * xi2 * nnu1 - nmu;
        nmu = nnu1;
        nnu1 = nnutemp;
      }
      nnu = nmu;
      npnu = nu * xi * nmu - nnu1;
      return kumi::make_tuple(jnu, jpnu, nnu, npnu);
    }
  } 
//   template<integral_scalar_value I, floating_real_scalar_value T>
//   EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I in, T x) noexcept
//   {
//     if( in == 0 )  return cyl_bessel_j0(x);
//     else  if( in == 1 )  return cyl_bessel_j1(x);
//     else return cyl_bessel_j(T(in), x);
//   }

//   template<integral_scalar_value I, floating_real_simd_value T>
//   EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I in, T x) noexcept
//   {
//     if( in == 0 )  return cyl_bessel_j0(x);
//     else  if( in == 1 )  return cyl_bessel_j1(x);
//     else return cyl_bessel_j(T(in), x);
//   }

//   template<integral_simd_value I, floating_real_simd_value T>
//   EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I in, T x) noexcept
//   requires compatible_values<T, I>
//   {
//     using elt_t =  element_type_t<T>;
//     return cyl_bessel_j(convert(in, as<elt_t>()), x);
//   }

//   template<floating_real_simd_value I, floating_real_scalar_value T>
//   EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), I in, T x) noexcept
//   requires compatible_values<T, I>
//   {
//     if constexpr( has_native_abi_v<T> )
//     {
//       return cyl_bessel_j(in, I(x));
//     }
//     else
//       return apply_over(cyl_bessel_j, in, x);
//   }

//   template<floating_real_value T>
//   EVE_FORCEINLINE auto cyl_bessel_j_(EVE_SUPPORTS(cpu_), T a0, T a1) noexcept
//   {
//     if constexpr( has_native_abi_v<T> )
//     {
//       auto sgn = if_else(is_ltz(a0), cospi(a0), one);
//       a0 = eve::abs(a0);
//       auto isneza1 =  is_nez(a1);
//       auto notdone = (!is_odd(a0) || is_nltz(a1)) && is_flint(a0) && isneza1;
//       T r = if_else(isneza1, nan(as(a1)), zero);
//       if (eve::any(notdone))
//       {
//         auto j0 = cyl_bessel_j0(a1);
//         auto br_0 =  [](auto j0) { return j0;};
//         notdone = next_interval(br_0, notdone, is_eqz(a0), r, j0);
//         if (eve::any(notdone))
//         {
//           auto j1 = cyl_bessel_j1(a1);
//           auto br_1 =  [](auto j1) { return j1;};
//           notdone = next_interval(br_1, notdone, a0 == 1, r, j1);
//           if (eve::any(notdone))
//           {
//             auto br_2 =  [](auto a1, auto j0,  auto j1)
//               {
//                 return fms(T(2)*j1, rec(a1), j0);
//               };
//             notdone = next_interval(br_2, notdone, a0 == 2, r, a1, j0, j1);
//             if (eve::any(notdone))
//             {
//               auto br_last = [](auto a0,  auto a1,  auto j0, auto j1)
//                 {
//                   std::int32_t k0 = 24;
//                   T pk = 2*(a0 + k0);
//                   auto ans = pk;
//                   auto xk = sqr(a1);
//                   do {
//                     pk  = pk - T(2);
//                     ans = fnma(xk, rec(ans), pk);
//                   }
//                   while( --k0 > 0 );
//                   /* backward recurrence */

//                   pk = T(1);
//                   /*pkm1 = 1.0/ans;*/
//                   T xinv = rec(a1);
//                   T pkm1 = ans * xinv;
//                   auto k = dec(a0);
//                   auto r = 2.0*k;
//                   auto test(true_(as(pk)));
//                   do{
//                     T pkm2 = (pkm1*r -  pk * a1) * xinv;
//                     pk   = if_else(test, pkm1, pk);
//                     pkm1 = if_else(test, pkm2, pkm1);
//                     r = r-T(2);
//                     k = dec(k);
//                     test = is_gtz(k);
//                   }
//                   while( eve::any(test) );
//                   return if_else(abs(pk) > pkm1, j1/pk, j0/pkm1);
//                 };
//               last_interval(br_last, notdone, r, a0, a1, j0, j1);
//             }
//           }
//         }
//       }
//       return sgn*r;
//     }
//     else
//       return apply_over(cyl_bessel_j, a0, a1);
//   }
}
