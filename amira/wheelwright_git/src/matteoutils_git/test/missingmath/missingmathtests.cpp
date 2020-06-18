/** @file missingmathtests.cpp
 * Tests for the functions in MissingMath.h and NumericTools.h
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
#include <limits>

#include "MissingMath.h"
#include "NumericTools.h"

#include "gtest/gtest.h"

TEST(MissingMathTest, binomWorks) {
    EXPECT_EQ(1ul, MatteoUtils::binom(1, 0));
    EXPECT_EQ(1ul, MatteoUtils::binom(1, 1));
    EXPECT_EQ(MatteoUtils::binom(5, 3), MatteoUtils::binom(5, 2));
    EXPECT_EQ(4950ul, MatteoUtils::binom(100, 2));
    EXPECT_EQ(1192052400ul, MatteoUtils::binom(100, 6));
}

TEST(MissingMathTest, logbinomWorks) {
    EXPECT_EQ(0, MatteoUtils::logbinom(1, 0));
    EXPECT_EQ(std::log(2), MatteoUtils::logbinom(2, 1));
    EXPECT_DOUBLE_EQ(MatteoUtils::logbinom(5, 3), MatteoUtils::logbinom(5, 2));
    EXPECT_DOUBLE_EQ(std::log(4950), MatteoUtils::logbinom(100, 2));
    EXPECT_DOUBLE_EQ(std::log(1192052400), MatteoUtils::logbinom(100, 6));
}

TEST(MissingMathTest, logcoshWorks) {
    EXPECT_EQ(0, MatteoUtils::logcosh(0));
    EXPECT_DOUBLE_EQ(19.306852819440056, MatteoUtils::logcosh(20));
    EXPECT_EQ(MatteoUtils::logcosh(20), MatteoUtils::logcosh(-20));
}


TEST(MissingMathTest, logsinhWorks) {
    EXPECT_DOUBLE_EQ(4.306807418479684, MatteoUtils::logsinh(5));
    using namespace std::complex_literals;
    EXPECT_DOUBLE_EQ(3.306530670702524, MatteoUtils::logsinh(4. + 3i).real());
    EXPECT_DOUBLE_EQ(2.9999062363415807, MatteoUtils::logsinh(4. + 3i).imag());
}

TEST(MissingMathTest, logsumlogWorks) {
    EXPECT_DOUBLE_EQ(std::log(7),
            MatteoUtils::logsumlog(std::log(4), std::log(3)));
}

TEST(MissingMathTest, logSumOfBinomsWorks) {
    EXPECT_DOUBLE_EQ(std::log(MatteoUtils::binom(100, 6)),
            MatteoUtils::log_sum_of_binoms(100, 6, 6));
    EXPECT_DOUBLE_EQ(std::log(1271427895ul),
            MatteoUtils::log_sum_of_binoms(100, 6));
    EXPECT_DOUBLE_EQ(std::log(1271427895ul - 100),
            MatteoUtils::log_sum_of_binoms(100, 6, 2));
}

TEST(MissingMathTest, sumOfBinomsWorks) {
    EXPECT_EQ(MatteoUtils::binom(100, 6),
            MatteoUtils::sum_of_binoms(100, 6, 6));
    EXPECT_EQ(1271427895ul, MatteoUtils::sum_of_binoms(100, 6));
}

TEST(MissingMathTest, logicallyEqualWorks) {
    EXPECT_TRUE(MatteoUtils::logically_equal(1, 1));
    EXPECT_FALSE(MatteoUtils::logically_equal(1, 2));
    EXPECT_TRUE(MatteoUtils::logically_equal(1.3, std::nexttoward(1.3, 1.31)));
    EXPECT_FALSE(MatteoUtils::logically_equal(1,
                1 + 2* std::numeric_limits<double>::epsilon()));
}
