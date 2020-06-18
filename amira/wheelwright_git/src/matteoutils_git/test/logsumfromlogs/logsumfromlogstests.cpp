/** @file logsumfromlogs.cpp
 * Tests for the class LogSumFromLogs
 *
 * @author Matteo Riondato
 * @date 2017 07 27
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

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "LogSumFromLogs.h"

#include "gtest/gtest.h"

TEST(LogSumFromLogsTest, LogSumFromLogsWorks) {
    MatteoUtils::LogSumFromLogs a {};
    EXPECT_DOUBLE_EQ(std::numeric_limits<double>::lowest(), a.value());
    MatteoUtils::LogSumFromLogs b {std::log(4)};
    EXPECT_DOUBLE_EQ(std::log(4), b.value());
    MatteoUtils::LogSumFromLogs c {std::log(5)};
    EXPECT_DOUBLE_EQ(std::log(9), (b + c).value());
    b += std::log(5);
    EXPECT_DOUBLE_EQ(std::log(9), b.value());
    c += b;
    EXPECT_DOUBLE_EQ(std::log(14), c.value());
    std::vector<MatteoUtils::LogSumFromLogs> v {b, c};
    #pragma omp declare reduction \
    (LogSumFromLogsPlus : MatteoUtils::LogSumFromLogs : omp_out += omp_in)
    #pragma omp parallel for reduction(LogSumFromLogsPlus : a)
    for (std::size_t i = 0; i < v.size(); ++i)
        a += v[i];
    EXPECT_DOUBLE_EQ(std::log(23), a.value());
}
