/** @file LogSumLog_impl.h
 * Implementation of the logsumlog function. It is not in MissingMath.h because
 * it is also used by the LogSumFromLog class, which is also used in
 * MissingMath.h by log_sum_of_binoms.
 *
 * @author Matteo Riondato
 * @date 2018 11 13
 *
 * @copyright
 *
 * Copyright 2018 Matteo Riondato <riondato@acm.org>
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

#ifndef MATTEOUTILS_LOGSUMLOGIMPL_H_
#define MATTEOUTILS_LOGSUMLOGIMPL_H_

#include <limits>
#include <utility>

namespace MatteoUtils::impl {

inline double logsumlog(double x, double y) {
    if (x <= std::numeric_limits<double>::lowest())
        return std::nexttoward(y, std::numeric_limits<double>::max());
    else if (y <= std::numeric_limits<double>::lowest())
        return std::nexttoward(x, std::numeric_limits<double>::max());
    else {
        if (x < y)
            std::swap(x, y);
        const double to_return {x + std::log1p(std::exp(y - x))};
        if (to_return < std::nexttoward(x, std::numeric_limits<double>::max()))
            return std::nexttoward(x, std::numeric_limits<double>::max());
        return to_return;
    }
}

} // namespace MatteoUtils::impl

#endif // MATTEOUTILS_LOGSUMLOGIMPL_H_
