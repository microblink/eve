//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

// Make MSVC compliant with macro we look for
#if defined(_MSC_VER)
#  if defined(EVE_ASSUME_SSE3)
#    define __SSE3__
#  endif

#  if defined(EVE_ASSUME_SSSE3)
#    define __SSSE3__
#  endif

#  if defined(EVE_ASSUME_SSE4_1)
#    define __SSE4_1__
#  endif

#  if defined(EVE_ASSUME_SSE4_2)
#    define __SSE4_2__
#  endif

#  if defined(EVE_ASSUME_XOP)
#    define __XOP__
#  endif

#  if defined(EVE_ASSUME_FMA4)
#    define __FMA4__
#  endif

#  if defined(EVE_ASSUME_FMA3)
#    define __FMA3__
#  endif
#endif

#include <eve/detail/spy.hpp>

// We successfully detected some native SIMD
#if defined(SPY_SIMD_IS_X86) && !defined(EVE_NO_SIMD)
#  define EVE_SUPPORTS_NATIVE_SIMD
#  define EVE_HW_X86
#  define EVE_INCLUDE_X86_HEADER

// Don't trigger AVX512 if we don't have at least Skylake support
#  if defined(SPY_SIMD_IS_X86_AVX512)
#   if !(   defined(__AVX512BW__) && defined(__AVX512CD__)  \
        &&  defined(__AVX512DQ__) && defined(__AVX512VL__)  \
        &&  defined(__AVX512DQ__) && defined(__AVX512VL__)  )
#    undef SPY_SIMD_IS_X86_AVX512
#    undef SPY_SIMD_DETECTED ::spy::detail::simd_version::avx512_
#    define EVE_INCOMPLETE_AVX512_SUPPORT
#    define SPY_SIMD_IS_X86_AVX2
#    define SPY_SIMD_DETECTED ::spy::detail::simd_version::avx2_
#   endif
#  endif

#ifdef __EMSCRIPTEN__
// SSE emulation on emscripten, doesn't support 256 and 512-bit AVX - add some stubs to silence compile errors
// in x86 headers that expect AVX stuff to exist

// Note: use #define instead of alias beacuse __m128 and friends are not visible yet

#define __m256d __m128d
#define __m256  __m128
#define __m256i __m128i

#define __m512d __m128d
#define __m512  __m128
#define __m512i __m128i

using __mmask8  = unsigned char;
using __mmask16 = unsigned short;
using __mmask32 = unsigned int;
using __mmask64 = unsigned long;

#define _mm256_movemask_pd _mm_movemask_pd
#define _mm256_movemask_ps _mm_movemask_ps

#define _mm512_load_si512  _mm_load_si128
#define _mm512_loadu_si512 _mm_loadu_si128

#define _mm256_load_si256  _mm_load_si128
#define _mm256_loadu_si256 _mm_loadu_si128

inline auto emptyFunc1Param( auto x )       noexcept { return x; }
inline auto emptyFunc2Param( auto x, auto ) noexcept { return x; }

#define _mm256_set_m128i       emptyFunc2Param
#define _mm256_shuffle_epi8    emptyFunc2Param
#define _mm256_zextsi128_si256 emptyFunc1Param
#define _mm256_permutevar8x32_epi32 emptyFunc2Param

#define _mm256_set_epi32 _mm_set_epi16
#define _mm256_set_epi16 _mm_set_epi8

#define _mm256_permute4x64_epi64 emptyFunc2Param
#define _mm256_permutevar8x32_ps emptyFunc2Param

#define _mm256_castps_pd _mm_castps_pd

#endif

#endif
