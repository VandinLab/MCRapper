/** @file MissingMath.h
 * Mathematical functions "missing" from the STL.
 *
 * @author Matteo Riondato
 * @date 2017 06 17
 *
 * @copyright
 *
 * Copyright 2017 Matteo Riondato <riondato@acm.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MATTEOUTILS_MISSINGMATH_H_
#define MATTEOUTILS_MISSINGMATH_H_

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#ifdef MATTEOUTILS_INVARIANTS
#include <stdexcept>
#endif

#include "LogSumFromLogs.h"

/** Main namespace for the library.
 *
 * All functions and classes belong to this namespace.
 *
 * Implementation details are in the non-documented  ``MatteoUtils::impl``
 * namespace.
 */
namespace MatteoUtils {

/** Precomputed natural logarithm of 2. */
constexpr static double log_2 { 0.69314718055994530941723212145818 };

/** Precomputed value for \f$\pi\f$. */
constexpr static double pi { 3.14159265358979323846 };

/** Computes the binomial coefficient \f$n \choose t\f$.
 *
 * @param n  the top value in the binomial coefficient.
 * @param t  the bottom value in the binomial coefficient.
 *
 * @returns \f$n \choose t\f$.
 */
inline std::size_t binom(const std::size_t n, const std::size_t t) {
    if (t == 0)
        return 1;
    if (t == 1)
        return n;
    if (2 * t > n)
        return binom(n, n - t);
#ifdef MATTEOUTILS_INVARIANTS
    if (t > n)
        throw std::invalid_argument("MatteoUtils::binom: t must be less than "
                "or equal to n.");
#endif
    std::size_t r {n - t + 1};
    for (std::size_t i {2}; i <= t; ++i) {
        r *= n - t + i;
        r /= i;
    }
    return r;
}

/** Computes the binomial coefficient \f$n \choose t\f$.
 *
 * @param n  the top value in the binomial coefficient.
 * @param t  the bottom value in the binomial coefficient.
 *
 * @returns \f$\ln n \choose t\f$.
 */
inline double logbinom(const std::size_t n, const std::size_t t) {
    if (t == 0)
        return 0;
    if (t == 1)
        return std::log(n);
    if (2 * t > n)
        return logbinom(n, n - t);
#ifdef MATTEOUTILS_INVARIANTS
    if (t > n)
        throw std::invalid_argument("MatteoUtils::logbinom: t must be less "
                "than or equal to n.");
#endif
    double r {std::log(n - t + 1)};
    for (std::size_t i {2}; i <= t; ++i)
        r += std::log(n - t + i) - std::log(i);
    return r;
}

/** Computes the natural logarithm of the hyperbolic cosine of @a x.
 *
 * @param x  the argument to the hyperbolic cosine.
 *
 * @returns \f$\ln\cosh(x)\f$.
 */
inline double logcosh(double x) {
    if (x == 0)
        return 0;
    else if (x < 0)
        x = -x;
    return x + std::log1p(std::exp(- 2 * x)) - log_2;
}

/** Computes the natural logarithm of the hyperbolic sine of a real number @a x.
 *
 * @param x the argument to the hyperbolic sine.
 *
 * @returns \f$\ln\sinh(x)\f$.
 */
inline double logsinh(const double x) {
#ifdef MATTEOUTILS_INVARIANTS
    if (x <= 0)
        throw std::invalid_argument("MatteoUtils::logsinh: the passed value is "
                "not positive");
#endif
    return x + std::log1p(- std::exp(- 2 * x)) - log_2;
}

/** Computes the natural logarithm of the hyperbolic sine of a complex number @a
 * x.
 *
 * @param x the argument to the hyperbolic sine.
 *
 * @returns \f$\ln\sinh(x)\f$, in the complex number sense.
 */
inline std::complex<double> logsinh(const std::complex<double> x) {
#ifdef MATTEOUTILS_INVARIANTS
    constexpr std::complex<double> zero;
    if (x == zero)
        throw std::invalid_argument("MatteoUtils::logsinh: the passed value is "
                "zero");
#endif
    return x + std::log(1.0 - std::exp(- 2.0 * x)) - log_2;
}

/** Computes the natural logarithm of the sum of the exponentials of the
 * arguments.
 *
 * @param x,y  the natural logarithms of the summands.
 *
 * @returns \f$\ln\left(e^x+e^y\right)\f$.
 */
inline double logsumlog(double x, double y) {
    return impl::logsumlog(x, y);
}

/** Computes the sum of binomial coefficients.
 *
 * @param n    the top value of the binomial coefficients.
 * @param max  the maximum bottom value of the binomial coefficients.
 * @param min  the minimum bottom value of the binomial coefficients.
 *
 * @returns \f$\sum_{i=min}^{max} {n \choose t}\f$.
 */
inline std::size_t sum_of_binoms(const std::size_t n, const std::size_t max,
        const std::size_t min = 1) {
#ifdef MATTEOUTILS_INVARIANTS
    if (max > n)
        throw std::out_of_range("MatteoUtils::sum_of_binoms: max must be at "
                "most n");
    if (min > max)
        throw std::out_of_range("MatteoUtils::sum_of_binoms: min must be at "
                "most max");
#endif
    std::size_t s {0};
    #pragma omp parallel for reduction(+: s)
    for (std::size_t i = min; i <= max; ++i)
        s += binom(n, i);
    return s;
}

/** Computes the logarithm of the sum of binomial coefficients.
 *
 * @param n    the top value of the binomial coefficients.
 * @param max  the maximum bottom value of the binomial coefficients.
 * @param min  the minimum bottom value of the binomial coefficients.
 *
 * @returns \f$\ln\sum_{i=min}^{max} {n \choose t}\f$.
 */

inline double log_sum_of_binoms(const std::size_t n, const std::size_t max,
        const std::size_t min = 1) {
#ifdef MATTEOUTILS_INVARIANTS
    if (max > n)
        throw std::out_of_range("MatteoUtils::log_sum_of_binoms: max must be "
                "at most n");
    if (min > max)
        throw std::out_of_range("MatteoUtils::log_sum_of_binoms: min must be "
                "at most max");
#endif
    // We use the recursive definition of the binomial coefficients.
    double coeff {0};
    for (std::size_t i {1}; i <= min; ++i)
        coeff += std::log(n - i + 1) - std::log(i);
    LogSumFromLogs sum {coeff};
    for (std::size_t i {min + 1}; i <= max ; ++i) {
        coeff += std::log(n - i + 1) - std::log(i);
        sum += coeff;
    }
    return sum.value();
}

} // namespace MatteoUtils

#endif // MATTEOUTILS_MISSINGMATH_H_
