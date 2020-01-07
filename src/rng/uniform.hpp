/*
Copyright 2016, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef rtt_r123_uniform_dot_hpp
#define rtt_r123_uniform_dot_hpp

// This file provides some simple tools that can be used to convert
// integers of various widths to floats and doubles with various
// characteristics.  It can be used to generate real-valued, uniformly
// distributed random variables from the random integers produced by
// the Random123 CBRNGs.
//
// There are three templated functions:
//
//   u01:  output is as dense as possible in (0,1}, never 0.0.  May
//     return 1.0 if and only if the number of output mantissa bits
//     is less than the width of the input.
//   uneg11:  output is as dense as possible in {-1,1}, never 0.0.  May
//     return 1.0 or -1.0 if and only if the number of output mantissa bits
//     is less than the width of the input.
//   u01fixedpt:  output is "fixed point", equispaced, open at both ends,
//      and is never 0.0, 0.5 nor 1.0.
//
// The behavior of u01 and uneg11 depend on the pre-processor symbol:
// R123_UNIFORM_FLOAT_STORE.  When #defined to a non-zero value, u01
// and uneg11 declare a volatile intermediate result, with the
// intention of forcing architectures that have "extra bits" in their
// floating point registers to more closely conform to IEEE
// arithmetic.  When compiled this way, u01 and uneg11 will be
// significantly slower, as they will incur a memory write and read on
// every call.  Without it, they may fail the "known answer test"
// implemented in ut_uniform_IEEEkat.cpp even though they perform
// perfectly reasonable int to float conversions.  We have used
// this option to get 32-bit x86 to produce the same results as
// 64-bit x86-64 code, but we do not recommend it for normal
// use.
//
// This file may not be as portable, and has not been tested as
// rigorously as the files in the library itself, i.e., those in
// ../include/Random123.  Nevertheless, we hope it is useful and we
// encourage developers to copy it and modify it for their own
// use.  We invite comments and improvements.

#include "rng/config.h"

#if defined(__GNUC__) && !defined(__clang__)
#if (DBS_GNUC_VERSION >= 70000)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wexpansion-to-defined"
#endif
#endif

#ifdef _MSC_FULL_VER
// conditional expression is constant
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <Random123/features/compilerfeatures.h>
#include <limits>
#if R123_USE_CXX11_TYPE_TRAITS
#include <type_traits>
#endif

namespace r123 {

#if R123_USE_CXX11_TYPE_TRAITS
using std::make_signed;
using std::make_unsigned;
#else
// Sigh... We could try to find another <type_traits>, e.g., from
// boost or TR1.  Or we can do it ourselves in the r123 namespace.
// It's not clear which will cause less headache...
template <typename T> struct make_signed {};
template <typename T> struct make_unsigned {};
#define R123_MK_SIGNED_UNSIGNED(ST, UT)                                        \
  template <> struct make_signed<ST> { typedef ST type; };                     \
  template <> struct make_signed<UT> { typedef ST type; };                     \
  template <> struct make_unsigned<ST> { typedef UT type; };                   \
  template <> struct make_unsigned<UT> { typedef UT type; }

R123_MK_SIGNED_UNSIGNED(int8_t, uint8_t);
R123_MK_SIGNED_UNSIGNED(int16_t, uint16_t);
R123_MK_SIGNED_UNSIGNED(int32_t, uint32_t);
R123_MK_SIGNED_UNSIGNED(int64_t, uint64_t);
#if R123_USE_GNU_UINT128
R123_MK_SIGNED_UNSIGNED(__int128_t, __uint128_t);
#endif
#undef R123_MK_SIGNED_UNSIGNED
#endif

#if defined(__CUDACC__) || defined(_LIBCPP_HAS_NO_CONSTEXPR)
// Amazing! cuda thinks numeric_limits::max() is a __host__ function, so
// we can't use it in a device function.
//
// The LIBCPP_HAS_NO_CONSTEXP test catches situations where the libc++
// library thinks that the compiler doesn't support constexpr, but we
// think it does.  As a consequence, the library declares
// numeric_limits::max without constexpr.  This workaround should only
// affect a narrow range of compiler/library pairings.
//
// In both cases, we find max() by computing ~(unsigned)0 right-shifted
// by is_signed.
template <typename T>
R123_CONSTEXPR R123_STATIC_INLINE R123_CUDA_DEVICE T maxTvalue() {
  typedef typename make_unsigned<T>::type uT;
  return (~uT(0)) >> std::numeric_limits<T>::is_signed;
}
#else
template <typename T> R123_CONSTEXPR R123_STATIC_INLINE T maxTvalue() {
  return std::numeric_limits<T>::max();
}
#endif

// u01: Input is a W-bit integer (signed or unsigned).  It is cast to
//   a W-bit unsigned integer, multiplied by Ftype(2^-W) and added to
//   Ftype(2^(-W-1)).  A good compiler should optimize it down to an
//   int-to-float conversion followed by a multiply and an add, which
//   might be fused, dependingon the architecture.
//
//  If the input is a uniformly distributed integer, then the
//  result is a uniformly distributed floating point number in [0, 1].
//  The result is never exactly 0.0.
//  The smallest value returned is 2^-W.
//  Let M be the number of mantissa bits in Ftype.
//  If W>M  then the largest value retured is 1.0.
//  If W<=M then the largest value returned is the largest Ftype less than 1.0.
template <typename Ftype, typename Itype>
R123_CUDA_DEVICE R123_STATIC_INLINE Ftype u01(Itype in) {
  typedef typename make_unsigned<Itype>::type Utype;
  R123_CONSTEXPR Ftype factor =
      Ftype(1.) / (static_cast<Ftype>(maxTvalue<Utype>()) + Ftype(1.));
  R123_CONSTEXPR Ftype halffactor = Ftype(0.5) * factor;
#ifdef R123_UNIFORM_FLOAT_STORE
  volatile Ftype x = Utype(in) * factor;
  return x + halffactor;
#else
  return Utype(in) * factor + halffactor;
#endif
}

// uneg11: Input is a W-bit integer (signed or unsigned).  It is cast
//    to a W-bit signed integer, multiplied by Ftype(2^-(W-1)) and
//    then added to Ftype(2^(-W-2)).  A good compiler should optimize
//    it down to an int-to-float conversion followed by a multiply and
//    an add, which might be fused, depending on the architecture.
//
//  If the input is a uniformly distributed integer, then the
//  output is a uniformly distributed floating point number in [-1, 1].
//  The result is never exactly 0.0.
//  The smallest absolute value returned is 2^-(W-1)
//  Let M be the number of mantissa bits in Ftype.
//  If W>M  then the largest value retured is 1.0 and the smallest is -1.0.
//  If W<=M then the largest value returned is the largest Ftype less than 1.0
//    and the smallest value returned is the smallest Ftype greater than -1.0.
template <typename Ftype, typename Itype>
R123_CUDA_DEVICE R123_STATIC_INLINE Ftype uneg11(Itype in) {
  typedef typename make_signed<Itype>::type Stype;
  R123_CONSTEXPR Ftype factor =
      Ftype(1.) / (Ftype(maxTvalue<Stype>()) + Ftype(1.));
  R123_CONSTEXPR Ftype halffactor = Ftype(0.5) * factor;
#ifdef R123_UNIFORM_FLOAT_STORE
  volatile Ftype x = Stype(in) * factor;
  return x + halffactor;
#else
  return Stype(in) * factor + halffactor;
#endif
}

// u01fixedpt:  Return a "fixed point" number in (0,1).  Let:
//   W = width of Itype, e.g., 32 or 64, regardless of signedness.
//   M = mantissa bits of Ftype, e.g., 24, 53 or 64
//   B = min(M, W)
// Then the 2^B-1 possible output values are:
//    2^-B*{1, 3, 5, ..., 2^B - 1}
// The smallest output is: 2^-B
// The largest output is:  1 - 2^-B
// The output is never exactly 0.0, nor 0.5, nor 1.0.
// The 2^(B-1) possible outputs:
//   - are equally likely,
//   - are uniformly spaced by 2^-(B-1),
//   - are balanced around 0.5
template <typename Ftype, typename Itype>
R123_CUDA_DEVICE R123_STATIC_INLINE Ftype u01fixedpt(Itype in) {
  typedef typename make_unsigned<Itype>::type Utype;
  R123_CONSTEXPR int excess =
      std::numeric_limits<Utype>::digits - std::numeric_limits<Ftype>::digits;
  if (excess >= 0) {

    // 2015-09-26 KT - Suppress warnings for the following expressions (see
    // https://rtt.lanl.gov/redmine/issues/416)
    //
    // Basically, GCC under BullseyeCoverage issues the following warning every
    // time this file is included:
    //
    // Counter_RNG.hh:124:65:   required from here
    // uniform.hpp:200:48: warning: second operand of conditional expression
    // has no effect [-Wunused-value]
    //         R123_CONSTEXPR int ex_nowarn = (excess>=0) ? excess : 0;
    //
    // Unfortunately, if this expression is simplified (see r7628) some
    // compilers will not compile the code because the RHS of the assignment
    // may contain values that are not known at comile time (not constexpr).
    // We don't want to spend to much time debugging this issue because the
    // code is essentially vendor owned (Random123).

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#endif

    R123_CONSTEXPR int ex_nowarn = (excess >= 0) ? excess : 0;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

    R123_CONSTEXPR Ftype factor =
        Ftype(1.) / (Ftype(1.) + ((maxTvalue<Utype>() >> ex_nowarn)));
    return (1 | (Utype(in) >> ex_nowarn)) * factor;
  } else {
    return u01<Ftype>(in);
  }
}

} // namespace r123

#ifdef _MSC_FULL_VER
// conditional expression is constant
#pragma warning(pop)
#endif

#ifdef __clang__
// Restore clang diagnostics to previous state.
#pragma clang diagnostic pop
#endif

#if defined(__GNUC__) && !defined(__clang__)
// Restore GCC diagnostics to previous state.
#pragma GCC diagnostic pop
#endif

#endif
