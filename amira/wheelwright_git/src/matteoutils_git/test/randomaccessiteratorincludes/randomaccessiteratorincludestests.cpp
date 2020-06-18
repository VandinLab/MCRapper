/** @file randomaccessiteratorincludestests.cpp
 * Tests for the function in RandomAccessIteratorIncludes.h
 *
 * @author Matteo Riondato
 * @date 2018 12 13
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

#include <vector>

#include "RandomAccessIteratorIncludes.h"

#include "gtest/gtest.h"

TEST(RandomAccessIteratorIncludesTest, includesWorks) {
    const std::vector<char> first {1, 2, 3, 4};
    EXPECT_TRUE(MatteoUtils::includes(first.cbegin(), first.cend(), first.cbegin(), first.cend()));
    const std::vector<char> sub {2, 4};
    EXPECT_TRUE(MatteoUtils::includes(first.cbegin(), first.cend(), sub.cbegin(), sub.cend()));
    const std::vector<char> empty;
    EXPECT_TRUE(MatteoUtils::includes(first.cbegin(), first.cend(), empty.cbegin(), empty.cend()));
    EXPECT_TRUE(MatteoUtils::includes(empty.cbegin(), empty.cend(), empty.cbegin(), empty.cend()));
    EXPECT_FALSE(MatteoUtils::includes(empty.cbegin(), empty.cend(), first.cbegin(), first.cend()));
    const std::vector<char> nsub1 {5, 6, 7};
    EXPECT_FALSE(MatteoUtils::includes(first.cbegin(), first.cend(), nsub1.cbegin(), nsub1.cend()));
    const std::vector<char> nsub2 {5, 6, 7, 8, 9};
    EXPECT_FALSE(MatteoUtils::includes(first.cbegin(), first.cend(), nsub1.cbegin(), nsub1.cend()));
}
