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
#include <eve/function/fms.hpp>
#include <eve/function/all.hpp>
#include <eve/function/bit_xor.hpp>
#include <eve/function/dec.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/is_infinite.hpp>
#include <eve/function/fma.hpp>
#include <eve/function/fms.hpp>
#include <eve/function/cyl_bessel_j1.hpp>
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
#include <eve/module/real/core/detail/generic/poleval.hpp>
#include <tuple>

namespace eve::detail
{

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
  EVE_FORCEINLINE V evaluate_rational(const T num, const U den, V z) noexcept
  {
    auto N =  num.size();
    auto eval_small = [&num, &den, N](auto z)
      {
        V s1(num[N-1]);
        V s2(den[N-1]);
        for(int i = (int)N - 2; i >= 0; --i)
        {
          s1 = fma(s1, z, num[i]);
          s2 = fma(s2, z, den[i]);
        }
        return s1/s2;
      };
    auto eval_large = [&num, &den, N](auto z)
      {
        z = rec(z);
        V s1(num[0]);
        V s2(den[0]);
        for(unsigned i = 1; i < N; ++i)
        {
          s1 = fma(s1, z, num[i]);
          s2 = fma(s2, z, den[i]);
        }
        return s1/s2;
       };
    auto test = z <= V(1);
    if(eve::all(test)) return eval_small(z);
    else if(eve::none(test)) return  eval_large(z);
    else return if_else(test,  eval_small(z),  eval_large(z));
  }


  template<floating_real_scalar_value T>
  EVE_FORCEINLINE T cyl_bessel_y1_(EVE_SUPPORTS(cpu_), T x) noexcept
  {
    std::cout << "cyl_bessel_y1_" << std::endl;
    auto Pi = eve::pi(as(x));
    if (x < 0)  return nan(as(x));
    if (x == 0) return minf(as(x));
    if (x == inf(as(x))) return zero(as(x));
    auto evaluate = [Pi](auto x, auto x1, auto x11, auto x12, auto P, auto Q)
      {
        T y = sqr(x);
        T z = 2 * log(x/x1) * cyl_bessel_j(T(1), x) / Pi;

//        T r = horner( y, P)/ horner( y, Q);
        T r = evaluate_rational(P, Q, y);
        T factor = (x + x1) * ((x - x11/256) - x12)/x;
        return fma(factor, r, z); //z + factor * r;
      };

    if (x <= T(4))                       // x in (0, 4]
    {
      constexpr std::array<T, 6> P1 = {
        4.0535726612579544093e+13
        ,  5.4708611716525426053e+12
        , -3.7595974497819597599e+11
        ,  7.2144548214502560419e+09
        , -5.9157479997408395984e+07
        ,  2.2157953222280260820e+05
        , -3.1714424660046133456e+02
      };
      constexpr std::array<T, 6> Q1 = {
        3.0737873921079286084e+14
        , 4.1272286200406461981e+12
        , 2.7800352738690585613e+10
        , 1.2250435122182963220e+08
        , 3.8136470753052572164e+05
        , 8.2079908168393867438e+02
        , 1.0
      };
      constexpr T
        x1  =   2.1971413260310170351e+00,
        x11 =   5.620e+02,
        x12 =   1.8288260310170351490e-03;
      return evaluate(x, x1, x11, x12, P1, Q1);
    }
    else if (x <= T(8))                  // x in (4, 8]
    {
      constexpr std::array<T, 7>  P2 = {
        ,  1.1514276357909013326e+19
        , -5.6808094574724204577e+18
        , -2.3638408497043134724e+16
        ,  4.0686275289804744814e+15
        , -5.9530713129741981618e+13
        ,  3.7453673962438488783e+11
        , -1.1957961912070617006e+09
        ,  1.9153806858264202986e+06
        , -1.2337180442012953128e+03
      };
      constexpr std::array<T, 7> Q2 = {
        , 5.3321844313316185697e+20
        , 5.6968198822857178911e+18
        , 3.0837179548112881950e+16
        , 1.1187010065856971027e+14
        , 3.0221766852960403645e+11
        , 6.3550318087088919566e+08
        , 1.0453748201934079734e+06
        , 1.2855164849321609336e+03
        , 1.0
      };
      constexpr T
        x2  =   5.4296810407941351328e+00,
        x21 =   1.3900e+03,
        x22 =   -6.4592058648672279948e-06;
      return evaluate(x, x2, x21, x22, P2, Q2);
    }
    else                                // x in (8, \infty)
    {
      constexpr std::array<T, 6> PC = {
         -4.4357578167941278571e+06
        , -9.9422465050776411957e+06
        , -6.6033732483649391093e+06
        , -1.5235293511811373833e+06
        , -1.0982405543459346727e+05
        , -1.6116166443246101165e+03
        , 0.0
      };
      constexpr std::array<T, 6> QC = {
         -4.4357578167941278568e+06
        , -9.9341243899345856590e+06
        , -6.5853394797230870728e+06
        , -1.5118095066341608816e+06
        , -1.0726385991103820119e+05
        , -1.4550094401904961825e+03
        , 1.0
      };
      constexpr std::array<T, 6> PS = {
          3.3220913409857223519e+04
         , 8.5145160675335701966e+04
         , 6.6178836581270835179e+04
         , 1.8494262873223866797e+04
         , 1.7063754290207680021e+03
         , 3.5265133846636032186e+01
         , 0.0
      };
      constexpr std::array<T, 6> QS = {
        7.0871281941028743574e+05
        , 1.8194580422439972989e+06
        , 1.4194606696037208929e+06
        , 4.0029443582266975117e+05
        , 3.7890229745772202641e+04
        , 8.6383677696049909675e+02
        , 1.0
      };
      T y = T(8)/x;
      T y2 = sqr(y);
      T rc = evaluate_rational(PC, QC, y2);//horner( y2, PC)/ horner( y2, QC);
      T rs = evaluate_rational(PS, QS, y2);//horner( y2, PS)/ horner( y2, QS);
      T factor = rsqrt(Pi*x);
      //
      // The following code is really just:
      //
      // T z = x - 0.75f * Pi
      // value = factor * (rc * sin(z) + y * rs * cos(z));
      //
      // But using the sin/cos addition formulae and constant values for
      // sin/cos of 3PI/4 which then cancel part of the "factor" term as they're all
      // 1 / sqrt(2):
      //
     auto [sx, cx] = sincos(x);
     return factor*fma(y, rs * (sx-cx), rc * (sx+cx));
    }
  }


  template<floating_real_simd_value T>
  EVE_FORCEINLINE T cyl_bessel_y1_(EVE_SUPPORTS(cpu_), T x) noexcept
  {
    if constexpr(has_native_abi_v<T>)
    {
      using v_t =  element_type_t<T>;
      auto Pi = eve::pi(as(x));
      T j1opi =  cyl_bessel_j(T(1), x) / Pi;

      auto evaluate = [j1opi](auto x, auto x1, auto x11, auto x12, auto P, auto Q)
        {
          T y = sqr(x);
          T z = 2 * log(x/x1) * j1opi;
//        T r = horner( y, P)/ horner( y, Q);
          T r = evaluate_rational(P, Q, y);
          T factor = (x + x1) * ((x - x11/256) - x12);
          return fma(factor, r, z);
        };

      auto br_4 = [evaluate](auto x){
        x =  if_else(x <= T(3), x, zero);
        constexpr std::array<v_t, 7> P1 = {
          4.0535726612579544093e+13
          ,  5.4708611716525426053e+12
          , -3.7595974497819597599e+11
          ,  7.2144548214502560419e+09
          , -5.9157479997408395984e+07
          ,  2.2157953222280260820e+05
          , -3.1714424660046133456e+02
        };
        constexpr std::array<v_t, 7> Q1 = {
          3.0737873921079286084e+14
          , 4.1272286200406461981e+12
          , 2.7800352738690585613e+10
          , 1.2250435122182963220e+08
          , 3.8136470753052572164e+05
          , 8.2079908168393867438e+02
          , 1.0
        };
        constexpr v_t
        x1  =   2.1971413260310170351e+00,
        x11 =   5.620e+02,
        x12 =   1.8288260310170351490e-03;
        return evaluate(x, x1, x11, x12, P1, Q1);
      };

      auto br_8 = [evaluate](auto x){
        constexpr std::array<v_t, 9>  P2 = {
          ,  1.1514276357909013326e+19
          , -5.6808094574724204577e+18
          , -2.3638408497043134724e+16
          ,  4.0686275289804744814e+15
          , -5.9530713129741981618e+13
          ,  3.7453673962438488783e+11
          , -1.1957961912070617006e+09
          ,  1.9153806858264202986e+06
          , -1.2337180442012953128e+03
        };
        constexpr std::array<v_t, 9> Q2 = {
          , 5.3321844313316185697e+20
          , 5.6968198822857178911e+18
          , 3.0837179548112881950e+16
          , 1.1187010065856971027e+14
          , 3.0221766852960403645e+11
          , 6.3550318087088919566e+08
          , 1.0453748201934079734e+06
          , 1.2855164849321609336e+03
          , 1.0
        };
        constexpr v_t
        x2  =   5.4296810407941351328e+00,
        x21 =   1.3900e+03,
        x22 =   -6.4592058648672279948e-06;
        return evaluate(x, x2, x21, x22, P2, Q2);
      };

      auto br_large = [Pi](auto x){
        constexpr std::array<v_t, 7> PC = {
          -4.4357578167941278571e+06
          , -9.9422465050776411957e+06
          , -6.6033732483649391093e+06
          , -1.5235293511811373833e+06
          , -1.0982405543459346727e+05
          , -1.6116166443246101165e+03
          , 0.0
        };
        constexpr std::array<v_t, 7> QC = {
          -4.4357578167941278568e+06
          , -9.9341243899345856590e+06
          , -6.5853394797230870728e+06
          , -1.5118095066341608816e+06
          , -1.0726385991103820119e+05
          , -1.4550094401904961825e+03
          , 1.0
        };
        constexpr std::array<v_t, 7> PS = {
          3.3220913409857223519e+04
          , 8.5145160675335701966e+04
          , 6.6178836581270835179e+04
          , 1.8494262873223866797e+04
          , 1.7063754290207680021e+03
          , 3.5265133846636032186e+01
          , 0.0
        };
        constexpr std::array<v_t, 7> QS = {
          7.0871281941028743574e+05
          , 1.8194580422439972989e+06
          , 1.4194606696037208929e+06
          , 4.0029443582266975117e+05
          , 3.7890229745772202641e+04
          , 8.6383677696049909675e+02
          , 1.0
        };
        T y = T(8)/x;
        T y2 = sqr(y);
        T rc = evaluate_rational(PC, QC, y2);//horner( y2, PC)/ horner( y2, QC);
        T rs = evaluate_rational(PS, QS, y2);//horner( y2, PS)/ horner( y2, QS);
        T factor = rsqrt(Pi*x);
        //
        // The following code is really just:
        //
        // T z = x - 0.75f * Pi
        // value = factor * (rc * sin(z) + y * rs * cos(z));
        //
        // But using the sin/cos addition formulae and constant values for
        // sin/cos of 3PI/4 which then cancel part of the "factor" term as they're all
        // 1 / sqrt(2):
        //
        auto [sx, cx] = sincos(x);
        return factor*fma(y, rs * (sx-cx), rc * (sx+cx));
      };

      auto r = nan(as(x));
      auto notdone = is_gez(x);

      if( eve::any(notdone) )
      {
        notdone = next_interval(br_3,  notdone, x <= T(3), r, x);
        if( eve::any(notdone) )
        {
          notdone = next_interval(br_5,  notdone, x <= T(5.5), r, x);
          if( eve::any(notdone) )
          {
            notdone = next_interval(br_8,  notdone, x <= T(8), r, x);
            if( eve::any(notdone) )
            {
              notdone = last_interval(br_large,  notdone, r, x);
            }
          }
        }
      }
      r = if_else (is_eqz(x), minf(as(x)), r);
      r = if_else (x == inf(as(x)), zero, r);
      return r;
    }
    else return apply_over(cyl_bessel_y0, x);
  }



 //  template<floating_real_value T>
//   EVE_FORCEINLINE T cyl_bessel_y1_(EVE_SUPPORTS(cpu_), T a0) noexcept
//   {
//     using elt_t =  element_type_t<T>;
//     if constexpr( has_native_abi_v<T> )
//     {
//       if constexpr(std::is_same_v<elt_t, float>)
//       {
//         auto branch1 =  [](auto x){
//           std::array<elt_t, 5> YP = {
//             8.061978323326852E-009f,
//             -9.496460629917016E-007f,
//             6.719543806674249E-005f,
//             -2.641785726447862E-003f,
//             4.202369946500099E-002f
//           };
//           constexpr elt_t Y01 = 4.66539330185668857532f;
//           const T z = sqr(x);
//           auto w = (z-Y01)*x*poleval( z, YP);
//           w = fma(twoopi(as(x)), fms(log(x), cyl_bessel_j1(x), rec(x)), w);
//           return if_else(is_eqz(x), minf(as(x)), w);
//         };

//         auto branch2 =  [](auto x){
//           std::array<elt_t, 8> MO1 = {
//             6.913942741265801E-002f,
//             -2.284801500053359E-001f,
//             3.138238455499697E-001f,
//             -2.102302420403875E-001f,
//             5.435364690523026E-003f,
//             1.493389585089498E-001f,
//             4.976029650847191E-006f,
//             7.978845453073848E-001f
//           };
//           std::array<elt_t, 8> PH1 = {
//             -4.497014141919556E+001f,
//             5.073465654089319E+001f,
//             -2.485774108720340E+001f,
//             7.222973196770240E+000f,
//             -1.544842782180211E+000f,
//             3.503787691653334E-001f,
//             -1.637986776941202E-001f,
//             3.749989509080821E-001f
//           };
//           auto q = rec(x);
//           auto w = sqrt(q);
//           auto p = w * poleval( q, MO1);
//           w = sqr(q);
//           auto xn = fms(q, poleval( w, PH1), 3*pio_4(as(x)));
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
//           std::array<elt_t, 6> YP = {
//             1.26320474790178026440E9,
//             -6.47355876379160291031E11,
//             1.14509511541823727583E14,
//             -8.12770255501325109621E15,
//             2.02439475713594898196E17,
//             -7.78877196265950026825E17,
//           };
//           std::array<elt_t, 8> YQ = {
//             /* 1.00000000000000000000E0,*/
//             5.94301592346128195359E2,
//             2.35564092943068577943E5,
//             7.34811944459721705660E7,
//             1.87601316108706159478E10,
//             3.88231277496238566008E12,
//             6.20557727146953693363E14,
//             6.87141087355300489866E16,
//             3.97270608116560655612E18,
//           };
//           const T z = sqr(x);
//           auto w = x*poleval( z, YP)/poleval1( z, YQ);
//           return if_else(is_eqz(x), minf(as(x))
//                         , w + twoopi(as(x)) * fms(log(x), cyl_bessel_j1(x), rec(x)));
//         };

//         auto branch2 =  [](auto x){
//           std::array<elt_t, 7> PP = {
//             7.62125616208173112003E-4,
//             7.31397056940917570436E-2,
//             1.12719608129684925192E0,
//             5.11207951146807644818E0,
//             8.42404590141772420927E0,
//             5.21451598682361504063E0,
//             1.00000000000000000254E0,
//           };
//           std::array<elt_t, 7>  PQ = {
//             5.71323128072548699714E-4,
//             6.88455908754495404082E-2,
//             1.10514232634061696926E0,
//             5.07386386128601488557E0,
//             8.39985554327604159757E0,
//             5.20982848682361821619E0,
//             9.99999999999999997461E-1,
//           };
//           std::array<elt_t, 8>  QP = {
//             5.10862594750176621635E-2,
//             4.98213872951233449420E0,
//             7.58238284132545283818E1,
//             3.66779609360150777800E2,
//             7.10856304998926107277E2,
//             5.97489612400613639965E2,
//             2.11688757100572135698E2,
//             2.52070205858023719784E1,
//           };
//           std::array<elt_t, 7> QQ = {
//             /*  1.00000000000000000000E0,*/
//             7.42373277035675149943E1,
//             1.05644886038262816351E3,
//             4.98641058337653607651E3,
//             9.56231892404756170795E3,
//             7.99704160447350683650E3,
//             2.82619278517639096600E3,
//             3.36093607810698293419E2,
//           };

//           auto  w = 5.0*rec(x);
//           auto  z = sqr(w);
//           auto  p = poleval( z, PP)/poleval( z, PQ );
//           auto  q = poleval( z, QP)/poleval1( z, QQ );
//           auto [s, c] = sincos(x);
//           p =  fms(w, q, p)*s - fma(w, q, p)*c;
//           p =  p*rsqrt(x*pi(as(x)));
//           return if_else(is_infinite(x), zero, p);
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
//       return apply_over(cyl_bessel_y1, a0);
//   }
}
