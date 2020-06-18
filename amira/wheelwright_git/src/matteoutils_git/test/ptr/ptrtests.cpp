/** @file ptrtests.cpp
 * Tests for the Ptr* classes.
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

#include <functional>
#include <memory>

#include "PtrEqualTo.h"
#include "PtrHash.h"
#include "PtrLess.h"

#include "gtest/gtest.h"

TEST(PtrTest, PtrEqualToWorks) {
    auto a = std::make_shared<int>(0);
    auto b = std::make_shared<int>(0);
    MatteoUtils::PtrEqualTo<int> t;
    EXPECT_TRUE(t(a, b));
    *b = 1;
    EXPECT_FALSE(t(a, b));
    EXPECT_TRUE(t(nullptr, nullptr));
    EXPECT_FALSE(t(a, nullptr));
    EXPECT_FALSE(t(nullptr, b));
}

TEST(PtrTest, PtrHashWorks) {
    auto a = std::make_shared<int>(0);
    MatteoUtils::PtrHash<int> t;
    EXPECT_EQ(std::hash<int>{}(0), t(a));
    EXPECT_NE(std::hash<int>{}(1), t(a));
}

TEST(PtrTest, PtrLessWorks) {
    auto a = std::make_shared<int>(0);
    auto b = std::make_shared<int>(1);
    MatteoUtils::PtrLess<int> t;
    EXPECT_TRUE(t(a, b));
    EXPECT_FALSE(t(b, a));
    EXPECT_FALSE(t(a, a));
    EXPECT_TRUE(t(nullptr, b));
    EXPECT_FALSE(t(nullptr, nullptr));
    EXPECT_FALSE(t(a, nullptr));
}
