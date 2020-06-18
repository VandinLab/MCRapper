/** @file LogSumFromLogs.h
 * The logarithm of a partial sum computed using the logarithms of the summands.
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

#ifndef MATTEOUTILS_LOGSUMFROMLOGS_H_
#define MATTEOUTILS_LOGSUMFROMLOGS_H_

#include <limits>
#include <stdexcept>

#include "LogSumLog_impl.h"

namespace MatteoUtils {

/** The logarithm of a partial sum computed using the logarithms of the
 * summands.
 *
 * This class is particularly useful to compute the sum in parallel, e.g.
 * using OpenMP using the reduction LogSumFromLogsPlus.
 */
class LogSumFromLogs {
    public:
        /** Default constructor with zero or one argument.
         *
         * @param init  initial value. Default: the most negative double.
         */
        LogSumFromLogs(const double init =
                std::numeric_limits<double>::lowest())
                : val {init} {
#ifdef MATTEOUTILS
            if (init <= 0)
                throw std::runtime_error(
                        "MatteoUtils::LogSumFromLogs::LogSumFromLogs: the "
                        "init value must be positive.");
#endif
        }

        /** Computes (in place) the sum with another LogSumFromLogs. */
        LogSumFromLogs& operator+=(const LogSumFromLogs &rhs) {
            this->val = impl::logsumlog(this->val, rhs.val);
            return *this;
        }

        /** Computes (in place) the sum with a double. */
        LogSumFromLogs& operator+=(const double rhs) {
#ifdef MATTEOUTILS
            if (rhs <= 0)
                throw std::runtime_error(
                        "MatteoUtils::LogSumFromLogs::operator+: the value to "
                        "add must be positive.");

#endif
            this->val = impl::logsumlog(this->val, rhs);
            return *this;
        }

        /** Set the value of the logarithm of the sum.
         *
         * @param v  the value to set it to.
         */
        void set(const double v) { val = v; }

        /** Returns the value of the logarithm of the sum. */
        double value() const { return val; }

        /** Computes the sum between two LogSumFromLogs objects.
         *
         * @param l,r  the two objects.
         *
         * @returns  a LogSumFromLogs object representing \f$\ln(e^l + e^r)\f$.
         */
        friend LogSumFromLogs operator+(LogSumFromLogs l,
                LogSumFromLogs &r) {
            l += r;
            return l;
        }
    private:
        double val;
}; // class LogSumFromLogs

/** An OpenMP reduction to compute partial logarithmic sum of an exponential
 * moments.
 */
#pragma omp declare reduction \
    (LogSumFromLogsPlus : LogSumFromLogs: omp_out += omp_in)

} // namespace MatteoUtils

#endif // MATTEOUTILS_LOGSUMFROMLOGS_H_
