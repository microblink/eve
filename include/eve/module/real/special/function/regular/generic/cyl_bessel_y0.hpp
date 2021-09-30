//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/detail/hz_device.hpp>
#include <eve/concept/value.hpp>
#include <eve/function/cyl_bessel_j.hpp>
#include <eve/function/fms.hpp>
#include <eve/function/all.hpp>
#include <eve/function/bit_xor.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/fma.hpp>
#include <eve/function/fms.hpp>
#include <eve/function/cyl_bessel_j0.hpp>
#include <eve/function/none.hpp>
#include <eve/function/sincos.hpp>
#include <eve/function/sin.hpp>
#include <eve/function/cos.hpp>
#include <eve/function/log.hpp>
#include <eve/function/rec.hpp>
#include <eve/function/rsqrt.hpp>
#include <eve/function/sqr.hpp>
#include <eve/function/sqrt.hpp>
#include <eve/constant/inf.hpp>
#include <eve/constant/invsqrt_2.hpp>
#include <eve/constant/pi.hpp>
#include <eve/constant/pio_4.hpp>
#include <type_traits>
#include <eve/function/horner.hpp>
#include <array>

namespace eve::detail
{

  template <class T, class U, class V>
  EVE_FORCEINLINE V evaluate_rational(const T num, const U denom, const V& z_) noexcept
  {
    auto N =  num.size();
    V z(z_);
    V s1, s2;
    if(z <= 1)
    {
      s1 = static_cast<V>(num[N-1]);
      s2 = static_cast<V>(denom[N-1]);
      for(int i = (int)N - 2; i >= 0; --i)
      {
        s1 *= z;
        s2 *= z;
        s1 += num[i];
        s2 += denom[i];
      }
    }
    else
    {
      z = rec(z);
      s1 = static_cast<V>(num[0]);
      s2 = static_cast<V>(denom[0]);
      for(unsigned i = 1; i < N; ++i)
      {
        s1 *= z;
        s2 *= z;
        s1 += num[i];
        s2 += denom[i];
      }
    }
    return s1 / s2;
  }


  template<floating_real_scalar_value T>
  EVE_FORCEINLINE T cyl_bessel_y0_(EVE_SUPPORTS(cpu_), T x) noexcept
  {
    std::cout << "cyl_bessel_y0_" << std::endl;
    constexpr std::array<T, 6> P1 = {
            1.0723538782003176831e+11
         , -8.3716255451260504098e+09
         ,  2.0422274357376619816e+08
         , -2.1287548474401797963e+06
         ,  1.0102532948020907590e+04
         , -1.8402381979244993524e+01
    };
    constexpr std::array<T, 6> Q1 = {
           5.8873865738997033405e+11
         , 8.1617187777290363573e+09
         , 5.5662956624278251596e+07
         , 2.3889393209447253406e+05
         , 6.6475986689240190091e+02
         , 1.0
    };
    constexpr std::array<T, 7>  P2 = {
          -2.2213976967566192242e+13
        , -5.5107435206722644429e+11
        ,  4.3600098638603061642e+10
        , -6.9590439394619619534e+08
        ,  4.6905288611678631510e+06
        , -1.4566865832663635920e+04
        ,  1.7427031242901594547e+01
    };
    constexpr std::array<T, 7> Q2 = {
           4.3386146580707264428e+14
         , 5.4266824419412347550e+12
         , 3.4015103849971240096e+10
         , 1.3960202770986831075e+08
         , 4.0669982352539552018e+05
         , 8.3030857612070288823e+02
         , 1.0
    };
    constexpr std::array<T, 8> P3 = {
         -8.0728726905150210443e+15
       ,  6.7016641869173237784e+14
       , -1.2829912364088687306e+11
       , -1.9363051266772083678e+11
       ,  2.1958827170518100757e+09
       , -1.0085539923498211426e+07
       ,  2.1363534169313901632e+04
       , -1.7439661319197499338e+01
    };
    constexpr std::array<T, 8> Q3 = {
           3.4563724628846457519e+17
         , 3.9272425569640309819e+15
         , 2.2598377924042897629e+13
         , 8.6926121104209825246e+10
         , 2.4727219475672302327e+08
         , 5.3924739209768057030e+05
         , 8.7903362168128450017e+02
         , 1.0
    };
    constexpr std::array<T, 6> PC = {
          2.2779090197304684302e+04
        , 4.1345386639580765797e+04
        , 2.1170523380864944322e+04
        , 3.4806486443249270347e+03
        , 1.5376201909008354296e+02
        , 8.8961548424210455236e-01
    };
    constexpr std::array<T, 6> QC = {
          2.2779090197304684318e+04
         , 4.1370412495510416640e+04
         , 2.1215350561880115730e+04
         , 3.5028735138235608207e+03
         , 1.5711159858080893649e+02
         , 1.0
    };
    constexpr std::array<T, 6> PS = {
          -8.9226600200800094098e+01
        , -1.8591953644342993800e+02
        , -1.1183429920482737611e+02
        , -2.2300261666214198472e+01
        , -1.2441026745835638459e+00
        , -8.8033303048680751817e-03
    };
    constexpr std::array<T, 6> QS = {
           5.7105024128512061905e+03
         , 1.1951131543434613647e+04
         , 7.2642780169211018836e+03
         , 1.4887231232283756582e+03
         , 9.0593769594993125859e+01
         , 1.0
    };
    constexpr T x1  =   8.9357696627916752158e-01,
      x2  =   3.9576784193148578684e+00,
      x3  =   7.0860510603017726976e+00,
      x11 =   2.280e+02,
      x12 =   2.9519662791675215849e-03,
      x21 =   1.0130e+03,
      x22 =   6.4716931485786837568e-04,
      x31 =   1.8140e+03,
      x32 =   1.1356030177269762362e-04;

 //    using namespace boost::math::tools;
//     using namespace boost::math::constants;

    auto Pi = eve::pi(as(x));
    if (x < 0)  return nan(as(x));
    if (x == 0) return minf(as(x));
    if (x == inf(as(x))) return zero(as(x));
    auto evaluate = [Pi](auto x, auto x1, auto x11, auto x12, auto P, auto Q)
      {
        T y = sqr(x);
        T z = 2 * log(x/x1) * cyl_bessel_j(T(0), x) / Pi;

//        T r = horner( y, P)/ horner( y, Q);
        T r = evaluate_rational(P, Q, y);
        T factor = (x + x1) * ((x - x11/256) - x12);
        return fma(factor, r, z); //z + factor * r;
      };

    if (x <= 3)                       // x in (0, 3]
    {
      return evaluate(x, x1, x11, x12, P1, Q1);
//       T y = sqr(x);
//       T z = 2 * log(x/x1) * cyl_bessel_j(T(0), x) / Pi;
//       T r = horner( y, P1)/ horner( y, Q1);
//  //      T r = evaluate_rational(P1, Q1, y);
//       T factor = (x + x1) * ((x - x11/256) - x12);
//       value = z + factor * r;
    }
    else if (x <= 5.5f)                  // x in (3, 5.5]
    {
      return evaluate(x, x2, x21, x22, P2, Q2);
//       T y = sqr(x);
//       T z = 2 * log(x/x2) * cyl_bessel_j(T(0), x) / Pi;
//       T r = evaluate_rational(P2, Q2, y);
//       T factor = (x + x2) * ((x - x21/256) - x22);
//      value = z + factor * r;
    }
    else if (x <= 8)                  // x in (5.5, 8]
    {
      return evaluate(x, x3, x31, x32, P3, Q3);
//       T y = sqr(x);
//       T z = 2 * log(x/x3) * cyl_bessel_j(T(0), x) / Pi;
//       T r = evaluate_rational(P3, Q3, y);
//       T factor = (x + x3) * ((x - x31/256) - x32);
//       value = z + factor * r;
    }
    else                                // x in (8, \infty)
    {
      std::cout << "here" << std::endl;
      T y = T(8)/x;
      T y2 = sqr(y);
      T rc = evaluate_rational(PC, QC, y2);//horner( y2, PC)/ horner( y2, QC);
      T rs = evaluate_rational(PS, QS, y2);//horner( y2, PS)/ horner( y2, QS);
      T factor =  T(5.641895835477562869480794515607725858e-01)/eve::sqrt(x); //rec(sqrt(Pi*x));
      //
      // The following code is really just:
      //
      // T z = x - 0.25f * Pi
      // value = factor * (rc * sin(z) + y * rs * cos(z));
      //
      // But using the sin/cos addition formulae and constant values for
      // sin/cos of PI/4 which then cancel part of the "factor" term as they're all
      // 1 / sqrt(2):
      //
//         T sx = sin(x);
//         T cx = cos(x);
     auto [sx, cx] = sincos(x);
        return factor*( rc * (sx - cx) +  y * rs * (cx + sx)) ;
    }
  }


  template<floating_real_simd_value T>
  EVE_FORCEINLINE T cyl_bessel_y0_(EVE_SUPPORTS(cpu_), T x) noexcept
  {
    return map(cyl_bessel_y0, x);
  }
}




//  *             cyl_bessel_y0()
//  * Bessel function of second kind, order zero  */
//  * Rational approximation coefficients YP[] are used for x < 6.5.
//  * The function computed is  cyl_bessel_y0(x)  -  2 ln(x) j0(x) / pi,
//  * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / pi
//  * = 0.073804295108687225 , EUL is Euler's constant.

//   template<floating_real_value T>
//   EVE_FORCEINLINE T cyl_bessel_y0_(EVE_SUPPORTS(cpu_), T a0) noexcept
//   {
//     using elt_t =  element_type_t<T>;
//     if constexpr( has_native_abi_v<T> )
//     {
//       if constexpr(std::is_same_v<elt_t, float>)
//       {
//         auto branch1 =  [](auto x){
//           std::array<elt_t, 5> YP = {
//             9.454583683980369E-008f,
//             -9.413212653797057E-006f,
//             5.344486707214273E-004f,
//             -1.584289289821316E-002f,
//             1.707584643733568E-001f
//           };
//           constexpr elt_t YZ1 =  0.43221455686510834878f;
//           const T z = sqr(x);
//           auto w = (z-YZ1) * poleval( z, YP);
//           return w + twoopi(as(x)) * log(x) *cyl_bessel_j0(x);
//         };

//         auto branch2 =  [](auto x){
//           std::array<elt_t, 8> MO = {
//             -6.838999669318810E-002f,
//             1.864949361379502E-001f,
//             -2.145007480346739E-001f,
//             1.197549369473540E-001f,
//             -3.560281861530129E-003f,
//             -4.969382655296620E-002f,
//             -3.355424622293709E-006f,
//             7.978845717621440E-001f
//           };
//           std::array<elt_t, 8> PH = {
//             3.242077816988247E+001f,
//             -3.630592630518434E+001f,
//             1.756221482109099E+001f,
//             -4.974978466280903E+000f,
//             1.001973420681837E+000f,
//             -1.939906941791308E-001f,
//             6.490598792654666E-002f,
//             -1.249992184872738E-001f
//           };
//           auto q = rec(x);
//           auto w = sqrt(q);
//           auto p = w * poleval( q, MO);
//           w = sqr(q);
//           auto xn = fms(q, poleval( w, PH), pio_4(as(x)));
//           return if_else(x == inf(as(x)), zero, p * sin(xn + x));
//         };

//         auto r = nan(as<T>()); //nan case treated here
//         auto notdone =  is_nltz(a0);
//         if(eve::any(notdone))
//         {
//           notdone = next_interval(branch1, notdone, a0 < T(2), r, a0);
//           if (eve::any(notdone))
//           {
//             last_interval(branch2, notdone, r, a0);
//           }
//         }
//         return r;
//       }
//       else // double
//       {
//         auto branch1 =  [](auto x){
//           std::array<elt_t, 8> YP = {
//             1.55924367855235737965E4,
//             -1.46639295903971606143E7,
//             5.43526477051876500413E9,
//             -9.82136065717911466409E11,
//             8.75906394395366999549E13,
//             -3.46628303384729719441E15,
//             4.42733268572569800351E16,
//             -1.84950800436986690637E16,
//           };
//           std::array<elt_t, 7> YQ = {
//             /* 1.00000000000000000000E0,*/
//             1.04128353664259848412E3,
//             6.26107330137134956842E5,
//             2.68919633393814121987E8,
//             8.64002487103935000337E10,
//             2.02979612750105546709E13,
//             3.17157752842975028269E15,
//             2.50596256172653059228E17,
//           };
//           const T z = sqr(x);
//           auto w = poleval( z, YP)/ poleval( z, YQ);
//           return w + twoopi(as(x)) * log(x) *cyl_bessel_j0(x);
//         };

//         auto branch2 =  [](auto x){
//           std::array<elt_t, 7> PP = {
//             7.96936729297347051624E-4,
//             8.28352392107440799803E-2,
//             1.23953371646414299388E0,
//             5.44725003058768775090E0,
//             8.74716500199817011941E0,
//             5.30324038235394892183E0,
//             9.99999999999999997821E-1,
//           };
//           std::array<elt_t, 7>  PQ = {
//             9.24408810558863637013E-4,
//             8.56288474354474431428E-2,
//             1.25352743901058953537E0,
//             5.47097740330417105182E0,
//             8.76190883237069594232E0,
//             5.30605288235394617618E0,
//             1.00000000000000000218E0,
//           };
//           std::array<elt_t, 8>  QP = {
//             -1.13663838898469149931E-2,
//             -1.28252718670509318512E0,
//             -1.95539544257735972385E1,
//             -9.32060152123768231369E1,
//             -1.77681167980488050595E2,
//             -1.47077505154951170175E2,
//             -5.14105326766599330220E1,
//             -6.05014350600728481186E0,
//           };
//           std::array<elt_t, 7> QQ = {
//             /*  1.00000000000000000000E0,*/
//             6.43178256118178023184E1,
//             8.56430025976980587198E2,
//             3.88240183605401609683E3,
//             7.24046774195652478189E3,
//             5.93072701187316984827E3,
//             2.06209331660327847417E3,
//             2.42005740240291393179E2,
//           };

//           auto  w = 5.0*rec(x);
//           auto  z = sqr(w);
//           auto  p = poleval( z, PP)/poleval( z, PQ );
//           auto  q = poleval( z, QP)/poleval1( z, QQ );
//           auto [s, c] = sincos(x);
//           p =  fma(w, q, p)*s +fms(w, q, p)*c;
//           return if_else(is_infinite(x), zero, p*rsqrt(x*pi(as(x))));
//         };

//         auto r = nan(as<T>()); //nan case treated here
//         auto notdone =  is_nltz(a0);
//         if(eve::any(notdone))
//         {
//           notdone = next_interval(branch1, notdone, a0 < T(5), r, a0);
//           if (eve::any(notdone))
//           {
//             last_interval(branch2, notdone, r, a0);
//           }
//         }
//         return r;
//       }
//     }
//     else
//       return apply_over(cyl_bessel_y0, a0);
//   }
// }
