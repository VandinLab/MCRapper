/** @file sizestets.cpp
 * Tests for the functions in VCDim.h
 *
 * @author Matteo Riondato
 * @date 2017 10 01
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
#include <cmath>
#include <list>
#include <vector>

#include "VCDim.h"

#include "gtest/gtest.h"

TEST(VCDimTest, VCDimWorks) {
    EXPECT_EQ(4922u, MatteoUtils::vcdim_size(0.05, 0.1, 10));
    EXPECT_DOUBLE_EQ(0.049995093406132778,
            MatteoUtils::vcdim_quality(0.1, 10, 4922u));

}
