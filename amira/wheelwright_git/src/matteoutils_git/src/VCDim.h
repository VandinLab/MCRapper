/** @file VCDim.h
 * Computation of sample size and maximum deviation bound using pseudodimension
 * or VC-dimension.
 *
 * @author Matteo Riondato
 * @date 2017 09 24
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

#ifndef MATTEOUTILS_VCDIM_H_
#define MATTEOUTILS_VCDIM_H_

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace MatteoUtils {

/** Compute a sample size to obtain an (@a eps, @a delta)-approximation to some
 * range space using an upper bound to the pseudodimension or the VC-dimension
 * of the range space.
 *
 * @param eps  the maximum allowed error.
 * @param d    the failure probability parameter (1 - confidence).
 * @param v    the upper bound to the pseudodimension or to the VC-dimension.
 */
std::size_t vcdim_size(const double eps, const double d,
        const std::size_t v) {
#ifdef MATTEOUTILS_INVARIANTS
    if (d <= 0 || d >=1)
        throw std::out_of_range("MatteoUtils::vcdim_quality(): The argument d "
                "must be in (0,1).");
    if (e <= 0 || e >=1)
        throw std::out_of_range("MatteoUtils::vcdim_size(): The argument e "
                "must be in (0,1).");
#endif
    // The 0.5 magic constant comes from a paper./
    // TODO MatteoUtils::vcdim_size: add the name of the paper.
    return static_cast<std::size_t>(std::ceil(0.5 * (v - std::log(d)) /
                std::pow(eps, 2)));
}

/** Compute a value that is, with probability at least 1 - @a d, an upper bound
 * to the maximum deviation of sample averages afrom their expectations, where
 * the averages are computed over a sample of size @a s and the VC-dimension or
 * pseudodimension of the problem are upper bounded by @a v.
 *
 * @param d    the failure probability parameter (1 - confidence).
 * @param v    the upper bound to the pseudodimension or to the VC-dimension.
 * @param s    the sample size.
 */
double vcdim_quality(const double d, const std::size_t v, const std::size_t s) {
    // The 2 magic constant has the same origin as the 0.5 in vcdim_size.
#ifdef MATTEOUTILS_INVARIANTS
    if (d <= 0 || d >=1)
        throw std::out_of_range("MatteoUtils::vcdim_quality(): The "
                "argument d must be in (0,1).");
#endif
    return std::sqrt((v - std::log(d)) / (2 * s));
}

} // namespace MatteoUtils

#endif // MATTEOUTILS_VCDIM_H_
