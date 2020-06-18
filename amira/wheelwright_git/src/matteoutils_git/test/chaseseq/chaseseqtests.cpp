/** @file chaseseqtests.cpp
 * Tests for the classes ChaseSeq and ChaseSeqIter
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

#include <cstddef>
#include <set>
#include <vector>

#include "ChaseSeq.h"

#include "gtest/gtest.h"

TEST(ChaseSeqTest, ChaseSeqWorks) {
    std::vector<std::size_t> v;
    constexpr std::size_t size {5};
    for (std::size_t i = 0; i < size; ++i)
        v.push_back(i);
    MatteoUtils::ChaseSeq<std::size_t> a(v.begin(), v.end(), 1,  1);
    EXPECT_EQ(size, a.size());
    std::set<std::set<std::size_t>> itms;
    for (auto s : a) {
        std::vector<std::size_t> v {s.first, s.second};
        EXPECT_EQ(1ul, v.size());
        itms.insert({v[0]});
    }
    EXPECT_EQ(size, itms.size());
    constexpr std::size_t size_choose_two {(size * (size - 1)) / 2};
    std::set<std::set<std::size_t>> pairs {{0,1}, {0,2}, {0,3}, {0,4}, {1,2},
        {1,3}, {1,4}, {2,3}, {2,4}, {3,4}};
    MatteoUtils::ChaseSeq<std::size_t> b(v.begin(), v.end(), 2, 2);
    EXPECT_EQ(size_choose_two, b.size());
    for (auto s : b) {
        std::set<std::size_t> v {s.first, s.second};
        EXPECT_EQ(2ul, v.size());
        EXPECT_TRUE(pairs.find(v) != pairs.end());
    }
    std::set<std::set<std::size_t>> others {{0,1,2}, {0,1,3}, {0,1,4}, {0,2,3},
        {0,2,4}, {0,3,4}, {1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}, {0,1,2,3},
        {0,1,2,4}, {0,1,3,4}, {0,2,3,4}, {1,2,3,4}, {0,1,2,3,4}};
    // TODO : MatteoUtils, chasesestests: We should test that the order in which
    // the sets are generated is correct.
    std::set<std::set<std::size_t>> all;
    all.insert(itms.begin(), itms.end());
    all.insert(pairs.begin(), pairs.end());
    all.insert(others.begin(), others.end());
    MatteoUtils::ChaseSeq<std::size_t> c(v.begin(), v.end(), size);
    EXPECT_EQ(all.size(), c.size());
    for (auto s : b) {
        std::set<std::size_t> v {s.first, s.second};
        EXPECT_TRUE(pairs.find(v) != pairs.end());
    }
}
