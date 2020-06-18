/** @file projectiontests.cpp
 * Tests for the classes Projection and ProjectionVecsIter
 *
 * @author Matteo Riondato
 * @date 2017 07 29
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
#include <iostream>
#include <set>
#include <utility>

#include "Projection.h"

#include "gtest/gtest.h"

TEST(ProjectionTest, ProjectionWorksWithVanilla) {
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p(7);
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    EXPECT_EQ(nullptr, p[0]);
    std::size_t i {0};
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test populating the Projection (simple case to then test drop_catchall
    // and clear. A more complex test is done later)
    p.begin_update(0, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(2ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(6ul, p.refs(p[7]));
    EXPECT_NE(p[7]->hash(), p[0]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // Test dims(d)
    p.dims(2);
    EXPECT_EQ(2ul, p.dims());
    // Test drop_catchall()
    p.drop_catchall();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(1ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(nullptr, p[7]);
    // Test clear()
    p.clear();
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    EXPECT_EQ(nullptr, p[0]);
    // Test starting to populate the Projection using the
    // begin_update()-update()-end_update() sequence.
    p.begin_update(0, 1);
    p.update(0, 1);
    p.update(1, 1);
    p.update(2, 1);
    p.update(3, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(0, p[0]->elltwosquared());
    EXPECT_EQ(2, p[3]->sum());
    EXPECT_EQ(1, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // update() of parameter that is the only parameter mapped to its vector
    p.begin_update(1, 1);
    p.update(3, 1);
    p.end_update();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0, p[7]->elltwosquared());
    // update() with new parameter causing a new vector to be created.
    p.begin_update(2, 2);
    p.update(4, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(4ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5, p[7]->elltwosquared());
    // Test update() with parameter which before the update is mapped to a
    // vector to which other parameters are mapped, and after the update is
    // mapped to a new vector.
    p.begin_update(4, 2);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(6ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[1], p[3]);
    EXPECT_NE(p[1], p[4]);
    EXPECT_NE(p[1], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[1]->hash(), p[3]->hash());
    EXPECT_NE(p[1]->hash(), p[4]->hash());
    EXPECT_NE(p[1]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[7]->elltwosquared());
    // Test centralized update dropping the cav vector when needed.
    p.begin_update(6, 1);
    p.update(5, 5);
    p.update(6, 5);
    p.end_update();
    EXPECT_EQ(7ul, p.dims());
    EXPECT_EQ(7ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(1ul, p.count(5));
    EXPECT_EQ(1ul, p.count(6));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[5]));
    EXPECT_EQ(2ul, p.refs(p[6]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[5]);
    EXPECT_NE(p[0], p[6]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[2], p[3]);
    EXPECT_NE(p[2], p[4]);
    EXPECT_NE(p[2], p[5]);
    EXPECT_NE(p[2], p[6]);
    EXPECT_NE(p[2], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[3], p[5]);
    EXPECT_NE(p[3], p[6]);
    EXPECT_NE(p[3], p[7]);
    EXPECT_NE(p[4], p[5]);
    EXPECT_NE(p[4], p[6]);
    EXPECT_NE(p[4], p[7]);
    EXPECT_EQ(p[5], p[6]);
    EXPECT_NE(p[5], p[7]);
    EXPECT_EQ(nullptr, p[7]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[5]->hash());
    EXPECT_NE(p[0]->hash(), p[6]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[2]->hash(), p[3]->hash());
    EXPECT_NE(p[2]->hash(), p[4]->hash());
    EXPECT_NE(p[2]->hash(), p[5]->hash());
    EXPECT_NE(p[2]->hash(), p[6]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[3]->hash(), p[5]->hash());
    EXPECT_NE(p[3]->hash(), p[6]->hash());
    EXPECT_NE(p[4]->hash(), p[5]->hash());
    EXPECT_NE(p[4]->hash(), p[6]->hash());
    EXPECT_EQ(p[5]->hash(), p[6]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0 + 1, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0 + 1, p[4]->elltwosquared());
    EXPECT_EQ(5, p[5]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 16, p[5]->elltwosquared());
    EXPECT_EQ(5, p[6]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 16, p[6]->elltwosquared());
    // Test merge()  -- 1
    // Merge two Projections such that:
    //  * one has the catch-all vector and the other doesn't;
    //  * there is a parameter that is mapped to the catch-all vector in one of
    //    them, and to a non-catch-all vector in the other.
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p1(2);
    p1.begin_update(0, 1);
    p1.update(0, 1);
    p1.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p2(2);
    p2.begin_update(1, 1);
    p2.update(0, 1);
    p2.update(1, 1);
    p2.end_update();
    p1.merge(p2);
    EXPECT_EQ(2ul, p1.dims());
    EXPECT_EQ(2ul, p1.par_size());
    EXPECT_EQ(2ul, p1.size());
    EXPECT_EQ(1ul, p1.count(0));
    EXPECT_EQ(1ul, p1.count(1));
    EXPECT_EQ(0ul, p1.count(2));
    i = 0;
    for (auto it {p1.cbegin()}; it != p1.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p1.refs(p1[0]));
    EXPECT_EQ(1ul, p1.refs(p1[1]));
    EXPECT_NE(p1[0], p1[1]);
    EXPECT_EQ(nullptr, p1[2]);
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(2, p1[0]->sum());
    EXPECT_DOUBLE_EQ(0.25, p1[0]->elltwosquared());
    EXPECT_EQ(1, p1[1]->sum());
    EXPECT_DOUBLE_EQ(0.25, p1[1]->elltwosquared());
    // Test merge() -- 2
    // Symmetric of the previous test.
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p3(2);
    p3.begin_update(1, 1);
    p3.update(0, 1);
    p3.update(1, 1);
    p3.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p4(2);
    p4.begin_update(0, 1);
    p4.update(0, 1);
    p4.end_update();
    p3.merge(p4);
    EXPECT_EQ(2ul, p3.dims());
    EXPECT_EQ(2ul, p3.par_size());
    EXPECT_EQ(2ul, p3.size());
    EXPECT_EQ(1ul, p3.count(0));
    EXPECT_EQ(1ul, p3.count(1));
    EXPECT_EQ(0ul, p3.count(2));
    i = 0;
    for (auto it {p3.cbegin()}; it != p3.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p3.refs(p3[0]));
    EXPECT_EQ(1ul, p3.refs(p3[1]));
    EXPECT_NE(p3[0], p3[1]);
    EXPECT_EQ(nullptr, p3[2]);
    EXPECT_NE(p3[0]->hash(), p3[1]->hash());
    EXPECT_EQ(2, p3[0]->sum());
    EXPECT_EQ(0.25, p3[0]->elltwosquared());
    EXPECT_EQ(1, p3[1]->sum());
    EXPECT_EQ(0.25, p3[1]->elltwosquared());
    // Test merge() -- 3
    // Merge two Projections both containing the catch-all vector but such that
    // the merged one does not contain the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p5(2);
    p5.begin_update(0, 1);
    p5.update(0, 1);
    p5.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p6(2);
    p6.begin_update(1, 1);
    p6.update(1, 1);
    p6.end_update();
    p5.merge(p6);
    EXPECT_EQ(2ul, p5.dims());
    EXPECT_EQ(2ul, p5.par_size());
    EXPECT_EQ(2ul, p5.size());
    EXPECT_EQ(1ul, p5.count(0));
    EXPECT_EQ(1ul, p5.count(1));
    EXPECT_EQ(0ul, p5.count(2));
    i = 0;
    for (auto it {p5.cbegin()}; it != p5.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p5.refs(p5[0]));
    EXPECT_EQ(1ul, p5.refs(p5[1]));
    EXPECT_NE(p5[0], p5[1]);
    EXPECT_EQ(nullptr, p5[2]);
    EXPECT_NE(p5[0]->hash(), p5[1]->hash());
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(0.5, p5[0]->elltwosquared());
    EXPECT_EQ(1, p5[1]->sum());
    EXPECT_EQ(0.5, p5[1]->elltwosquared());
    // Test merge() -- 4
    // Merge two Projections both containing the catch-all vector and such that
    // the merged one also contains the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p7(2);
    p7.begin_update(0, 1);
    p7.update(0, 1);
    p7.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p8(2);
    p8.begin_update(1, 1);
    p8.update(0, 1);
    p8.end_update();
    p7.merge(p8);
    EXPECT_EQ(2ul, p7.dims());
    EXPECT_EQ(1ul, p7.par_size());
    EXPECT_EQ(2ul, p7.size());
    EXPECT_EQ(1ul, p7.count(0));
    EXPECT_EQ(0ul, p7.count(1));
    EXPECT_EQ(0ul, p7.count(2));
    i = 0;
    for (auto it {p7.cbegin()}; it != p7.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p7.refs(p7[0]));
    EXPECT_EQ(1ul, p7.refs(p7[2]));
    EXPECT_NE(p7[0], p7[2]);
    EXPECT_EQ(2, p7[0]->sum());
    EXPECT_EQ(0.5, p7[0]->elltwosquared());
    EXPECT_EQ(0, p7[2]->sum());
    EXPECT_EQ(0.5, p7[2]->elltwosquared());
    // Test merge()  -- 5
    // Merge two Projections such that a pair of parameters is mapped to the
    // same vector in one Projection and to two different vectors in the other.
    // Vice versa for another pair of parameters.
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p9(6);
    p9.begin_update(0, 1);
    p9.update(0, 1);
    p9.update(1, 1);
    p9.update(2, 2);
    p9.update(3, 3);
    p9.update(4, 3);
    p9.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p10(6);
    p10.begin_update(1, 1);
    p10.update(0, 1);
    p10.update(1, 2);
    p10.update(2, 3);
    p10.update(3, 3);
    p10.end_update();
    p9.merge(p10);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge()  -- 6
    // Test the special case that the oth Projection is empty.
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> e(6);
    p9.merge(e);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge() -- 7
    // Test the special case that this Projection is empty.
    e.merge(p9);
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy constructor;
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p11 {e};
    EXPECT_EQ(2ul, p11.dims());
    EXPECT_EQ(5ul, p11.par_size());
    EXPECT_EQ(6ul, p11.size());
    EXPECT_EQ(1ul, p11.count(0));
    EXPECT_EQ(1ul, p11.count(1));
    EXPECT_EQ(1ul, p11.count(2));
    EXPECT_EQ(1ul, p11.count(3));
    EXPECT_EQ(1ul, p11.count(4));
    EXPECT_EQ(0ul, p11.count(5));
    EXPECT_EQ(0ul, p11.count(6));
    i = 0;
    for (auto it {p11.cbegin()}; it != p11.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p11.refs(p11[0]));
    EXPECT_EQ(1ul, p11.refs(p11[1]));
    EXPECT_EQ(1ul, p11.refs(p11[2]));
    EXPECT_EQ(1ul, p11.refs(p11[3]));
    EXPECT_EQ(1ul, p11.refs(p11[4]));
    EXPECT_EQ(1ul, p11.refs(p11[6]));
    EXPECT_EQ(2, p11[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p11[0]->elltwosquared());
    EXPECT_EQ(3, p11[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p11[1]->elltwosquared());
    EXPECT_EQ(5, p11[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p11[2]->elltwosquared());
    EXPECT_EQ(6, p11[3]->sum());
    EXPECT_EQ(4.5, p11[3]->elltwosquared());
    EXPECT_EQ(3, p11[4]->sum());
    EXPECT_EQ(4.5, p11[4]->elltwosquared());
    EXPECT_EQ(0, p11[6]->sum());
    EXPECT_EQ(4.5, p11[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.par_cbegin()}; it != e.par_cend(); ++it)
        ++i;
    EXPECT_EQ(e.par_size(), i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy assignment
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA>  p12(6);
    p12 = e;
    EXPECT_EQ(2ul, p12.dims());
    EXPECT_EQ(5ul, p12.par_size());
    EXPECT_EQ(6ul, p12.size());
    EXPECT_EQ(1ul, p12.count(0));
    EXPECT_EQ(1ul, p12.count(1));
    EXPECT_EQ(1ul, p12.count(2));
    EXPECT_EQ(1ul, p12.count(3));
    EXPECT_EQ(1ul, p12.count(4));
    EXPECT_EQ(0ul, p12.count(5));
    EXPECT_EQ(0ul, p12.count(6));
    i = 0;
    for (auto it {p12.cbegin()}; it != p12.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p12.refs(p12[0]));
    EXPECT_EQ(1ul, p12.refs(p12[1]));
    EXPECT_EQ(1ul, p12.refs(p12[2]));
    EXPECT_EQ(1ul, p12.refs(p12[3]));
    EXPECT_EQ(1ul, p12.refs(p12[4]));
    EXPECT_EQ(1ul, p12.refs(p12[6]));
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p12[0]->elltwosquared());
    EXPECT_EQ(3, p12[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p12[1]->elltwosquared());
    EXPECT_EQ(5, p12[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p12[2]->elltwosquared());
    EXPECT_EQ(6, p12[3]->sum());
    EXPECT_EQ(4.5, p12[3]->elltwosquared());
    EXPECT_EQ(3, p12[4]->sum());
    EXPECT_EQ(4.5, p12[4]->elltwosquared());
    EXPECT_EQ(0, p12[6]->sum());
    EXPECT_EQ(4.5, p12[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test move constructor
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p13 {std::move(e)};
    EXPECT_EQ(2ul, p13.dims());
    EXPECT_EQ(5ul, p13.par_size());
    EXPECT_EQ(6ul, p13.size());
    EXPECT_EQ(1ul, p13.count(0));
    EXPECT_EQ(1ul, p13.count(1));
    EXPECT_EQ(1ul, p13.count(2));
    EXPECT_EQ(1ul, p13.count(3));
    EXPECT_EQ(1ul, p13.count(4));
    EXPECT_EQ(0ul, p13.count(5));
    EXPECT_EQ(0ul, p13.count(6));
    i = 0;
    for (auto it {p13.cbegin()}; it != p13.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p13.refs(p13[0]));
    EXPECT_EQ(1ul, p13.refs(p13[1]));
    EXPECT_EQ(1ul, p13.refs(p13[2]));
    EXPECT_EQ(1ul, p13.refs(p13[3]));
    EXPECT_EQ(1ul, p13.refs(p13[4]));
    EXPECT_EQ(1ul, p13.refs(p13[6]));
    EXPECT_EQ(2, p13[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p13[0]->elltwosquared());
    EXPECT_EQ(3, p13[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p13[1]->elltwosquared());
    EXPECT_EQ(5, p13[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p13[2]->elltwosquared());
    EXPECT_EQ(6, p13[3]->sum());
    EXPECT_EQ(4.5, p13[3]->elltwosquared());
    EXPECT_EQ(3, p13[4]->sum());
    EXPECT_EQ(4.5, p13[4]->elltwosquared());
    EXPECT_EQ(0, p13[6]->sum());
    EXPECT_EQ(4.5, p13[6]->elltwosquared());
    // Test move assignment
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA> p14(6);
    p14 = std::move(p12);
    EXPECT_EQ(2ul, p14.dims());
    EXPECT_EQ(5ul, p14.par_size());
    EXPECT_EQ(6ul, p14.size());
    EXPECT_EQ(1ul, p14.count(0));
    EXPECT_EQ(1ul, p14.count(1));
    EXPECT_EQ(1ul, p14.count(2));
    EXPECT_EQ(1ul, p14.count(3));
    EXPECT_EQ(1ul, p14.count(4));
    EXPECT_EQ(0ul, p14.count(5));
    EXPECT_EQ(0ul, p14.count(6));
    i = 0;
    for (auto it {p14.cbegin()}; it != p14.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p14.refs(p14[0]));
    EXPECT_EQ(1ul, p14.refs(p14[1]));
    EXPECT_EQ(1ul, p14.refs(p14[2]));
    EXPECT_EQ(1ul, p14.refs(p14[3]));
    EXPECT_EQ(1ul, p14.refs(p14[4]));
    EXPECT_EQ(1ul, p14.refs(p14[6]));
    EXPECT_EQ(2, p14[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p14[0]->elltwosquared());
    EXPECT_EQ(3, p14[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p14[1]->elltwosquared());
    EXPECT_EQ(5, p14[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p14[2]->elltwosquared());
    EXPECT_EQ(6, p14[3]->sum());
    EXPECT_EQ(4.5, p14[3]->elltwosquared());
    EXPECT_EQ(3, p14[4]->sum());
    EXPECT_EQ(4.5, p14[4]->elltwosquared());
    EXPECT_EQ(0, p14[6]->sum());
    EXPECT_EQ(4.5, p14[6]->elltwosquared());
}


TEST(ProjectionTest, ProjectionWorksWithDiscreteAndVanilla) {
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    std::set<int> v {1, 2, 3, 4, 5};
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p(7, v.begin(), v.end());
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    EXPECT_EQ(nullptr, p[0]);
    std::size_t i {0};
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test populating the Projection (simple case to then test drop_catchall
    // and clear. A more complex test is done later)
    p.begin_update(0, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(2ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(6ul, p.refs(p[7]));
    EXPECT_NE(p[7]->hash(), p[0]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // Test dims(d)
    p.dims(2);
    EXPECT_EQ(2ul, p.dims());
    // Test drop_catchall()
    p.drop_catchall();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(1ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(nullptr, p[7]);
    // Test clear()
    p.clear();
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    EXPECT_EQ(nullptr, p[0]);
    // Test starting to populate the Projection using the
    // begin_update()-update()-end_update() sequence.
    p.begin_update(0, 1);
    p.update(0, 1);
    p.update(1, 1);
    p.update(2, 1);
    p.update(3, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(0, p[0]->elltwosquared());
    EXPECT_EQ(2, p[3]->sum());
    EXPECT_EQ(1, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // update() of parameter that is the only parameter mapped to its vector
    p.begin_update(1, 1);
    p.update(3, 1);
    p.end_update();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0, p[7]->elltwosquared());
    // update() with new parameter causing a new vector to be created.
    p.begin_update(2, 2);
    p.update(4, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(4ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5, p[7]->elltwosquared());
    // Test update() with parameter which before the update is mapped to a
    // vector to which other parameters are mapped, and after the update is
    // mapped to a new vector.
    p.begin_update(4, 2);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(6ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[1], p[3]);
    EXPECT_NE(p[1], p[4]);
    EXPECT_NE(p[1], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[1]->hash(), p[3]->hash());
    EXPECT_NE(p[1]->hash(), p[4]->hash());
    EXPECT_NE(p[1]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[7]->elltwosquared());
    // Test centralized update dropping the cav vector when needed.
    p.begin_update(6, 1);
    p.update(5, 5);
    p.update(6, 5);
    p.end_update();
    EXPECT_EQ(7ul, p.dims());
    EXPECT_EQ(7ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(1ul, p.count(5));
    EXPECT_EQ(1ul, p.count(6));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[5]));
    EXPECT_EQ(2ul, p.refs(p[6]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[5]);
    EXPECT_NE(p[0], p[6]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[2], p[3]);
    EXPECT_NE(p[2], p[4]);
    EXPECT_NE(p[2], p[5]);
    EXPECT_NE(p[2], p[6]);
    EXPECT_NE(p[2], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[3], p[5]);
    EXPECT_NE(p[3], p[6]);
    EXPECT_NE(p[3], p[7]);
    EXPECT_NE(p[4], p[5]);
    EXPECT_NE(p[4], p[6]);
    EXPECT_NE(p[4], p[7]);
    EXPECT_EQ(p[5], p[6]);
    EXPECT_NE(p[5], p[7]);
    EXPECT_EQ(nullptr, p[7]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[5]->hash());
    EXPECT_NE(p[0]->hash(), p[6]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[2]->hash(), p[3]->hash());
    EXPECT_NE(p[2]->hash(), p[4]->hash());
    EXPECT_NE(p[2]->hash(), p[5]->hash());
    EXPECT_NE(p[2]->hash(), p[6]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[3]->hash(), p[5]->hash());
    EXPECT_NE(p[3]->hash(), p[6]->hash());
    EXPECT_NE(p[4]->hash(), p[5]->hash());
    EXPECT_NE(p[4]->hash(), p[6]->hash());
    EXPECT_EQ(p[5]->hash(), p[6]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0 + 1, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0 + 1, p[4]->elltwosquared());
    EXPECT_EQ(5, p[5]->sum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 16, p[5]->elltwosquared());
    // Test merge()  -- 1
    // Merge two Projections such that:
    //  * one has the catch-all vector and the other doesn't;
    //  * there is a parameter that is mapped to the catch-all vector in one of
    //    them, and to a non-catch-all vector in the other.
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p1(2, v.begin(), v.end());
    p1.begin_update(0, 1);
    p1.update(0, 1);
    p1.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p2(2, v.begin(), v.end());
    p2.begin_update(1, 1);
    p2.update(0, 1);
    p2.update(1, 1);
    p2.end_update();
    p1.merge(p2);
    EXPECT_EQ(2ul, p1.dims());
    EXPECT_EQ(2ul, p1.par_size());
    EXPECT_EQ(2ul, p1.size());
    EXPECT_EQ(1ul, p1.count(0));
    EXPECT_EQ(1ul, p1.count(1));
    EXPECT_EQ(0ul, p1.count(2));
    i = 0;
    for (auto it {p1.cbegin()}; it != p1.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p1.refs(p1[0]));
    EXPECT_EQ(1ul, p1.refs(p1[1]));
    EXPECT_NE(p1[0], p1[1]);
    EXPECT_EQ(nullptr, p1[2]);
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(2, p1[0]->sum());
    EXPECT_DOUBLE_EQ(0.25, p1[0]->elltwosquared());
    EXPECT_EQ(1, p1[1]->sum());
    EXPECT_DOUBLE_EQ(0.25, p1[1]->elltwosquared());
    // Test merge() -- 2
    // Symmetric of the previous test.
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p3(2, v.begin(), v.end());
    p3.begin_update(1, 1);
    p3.update(0, 1);
    p3.update(1, 1);
    p3.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p4(2, v.begin(), v.end());
    p4.begin_update(0, 1);
    p4.update(0, 1);
    p4.end_update();
    p3.merge(p4);
    EXPECT_EQ(2ul, p3.dims());
    EXPECT_EQ(2ul, p3.par_size());
    EXPECT_EQ(2ul, p3.size());
    EXPECT_EQ(1ul, p3.count(0));
    EXPECT_EQ(1ul, p3.count(1));
    EXPECT_EQ(0ul, p3.count(2));
    i = 0;
    for (auto it {p3.cbegin()}; it != p3.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p3.refs(p3[0]));
    EXPECT_EQ(1ul, p3.refs(p3[1]));
    EXPECT_NE(p3[0], p3[1]);
    EXPECT_EQ(nullptr, p3[2]);
    EXPECT_NE(p3[0]->hash(), p3[1]->hash());
    EXPECT_EQ(2, p3[0]->sum());
    EXPECT_EQ(0.25, p3[0]->elltwosquared());
    EXPECT_EQ(1, p3[1]->sum());
    EXPECT_EQ(0.25, p3[1]->elltwosquared());
    // Test merge() -- 3
    // Merge two Projections both containing the catch-all vector but such that
    // the merged one does not contain the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p5(2, v.begin(), v.end());
    p5.begin_update(0, 1);
    p5.update(0, 1);
    p5.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p6(2, v.begin(), v.end());
    p6.begin_update(1, 1);
    p6.update(1, 1);
    p6.end_update();
    p5.merge(p6);
    EXPECT_EQ(2ul, p5.dims());
    EXPECT_EQ(2ul, p5.par_size());
    EXPECT_EQ(2ul, p5.size());
    EXPECT_EQ(1ul, p5.count(0));
    EXPECT_EQ(1ul, p5.count(1));
    EXPECT_EQ(0ul, p5.count(2));
    i = 0;
    for (auto it {p5.cbegin()}; it != p5.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p5.refs(p5[0]));
    EXPECT_EQ(1ul, p5.refs(p5[1]));
    EXPECT_NE(p5[0], p5[1]);
    EXPECT_EQ(nullptr, p5[2]);
    EXPECT_NE(p5[0]->hash(), p5[1]->hash());
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(0.5, p5[0]->elltwosquared());
    EXPECT_EQ(1, p5[1]->sum());
    EXPECT_EQ(0.5, p5[1]->elltwosquared());
    // Test merge() -- 4
    // Merge two Projections both containing the catch-all vector and such that
    // the merged one also contains the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p7(2, v.begin(), v.end());
    p7.begin_update(0, 1);
    p7.update(0, 1);
    p7.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p8(2, v.begin(), v.end());
    p8.begin_update(1, 1);
    p8.update(0, 1);
    p8.end_update();
    p7.merge(p8);
    EXPECT_EQ(2ul, p7.dims());
    EXPECT_EQ(1ul, p7.par_size());
    EXPECT_EQ(2ul, p7.size());
    EXPECT_EQ(1ul, p7.count(0));
    EXPECT_EQ(0ul, p7.count(1));
    EXPECT_EQ(0ul, p7.count(2));
    i = 0;
    for (auto it {p7.cbegin()}; it != p7.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p7.refs(p7[0]));
    EXPECT_EQ(1ul, p7.refs(p7[2]));
    EXPECT_NE(p7[0], p7[2]);
    EXPECT_EQ(2, p7[0]->sum());
    EXPECT_EQ(0.5, p7[0]->elltwosquared());
    EXPECT_EQ(0, p7[2]->sum());
    EXPECT_EQ(0.5, p7[2]->elltwosquared());
    // Test merge()  -- 5
    // Merge two Projections such that a pair of parameters is mapped to the
    // same vector in one Projection and to two different vectors in the other.
    // Vice versa for another pair of parameters.
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p9(6, v.begin(), v.end());
    p9.begin_update(0, 1);
    p9.update(0, 1);
    p9.update(1, 1);
    p9.update(2, 2);
    p9.update(3, 3);
    p9.update(4, 3);
    p9.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p10(6, v.begin(), v.end());
    p10.begin_update(1, 1);
    p10.update(0, 1);
    p10.update(1, 2);
    p10.update(2, 3);
    p10.update(3, 3);
    p10.end_update();
    p9.merge(p10);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge()  -- 6
    // Test the special case that the oth Projection is empty.
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> e(6, v.begin(), v.end());
    p9.merge(e);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge() -- 7
    // Test the special case that this Projection is empty.
    e.merge(p9);
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy constructor;
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p11 {e};
    EXPECT_EQ(2ul, p11.dims());
    EXPECT_EQ(5ul, p11.par_size());
    EXPECT_EQ(6ul, p11.size());
    EXPECT_EQ(1ul, p11.count(0));
    EXPECT_EQ(1ul, p11.count(1));
    EXPECT_EQ(1ul, p11.count(2));
    EXPECT_EQ(1ul, p11.count(3));
    EXPECT_EQ(1ul, p11.count(4));
    EXPECT_EQ(0ul, p11.count(5));
    EXPECT_EQ(0ul, p11.count(6));
    i = 0;
    for (auto it {p11.cbegin()}; it != p11.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p11.refs(p11[0]));
    EXPECT_EQ(1ul, p11.refs(p11[1]));
    EXPECT_EQ(1ul, p11.refs(p11[2]));
    EXPECT_EQ(1ul, p11.refs(p11[3]));
    EXPECT_EQ(1ul, p11.refs(p11[4]));
    EXPECT_EQ(1ul, p11.refs(p11[6]));
    EXPECT_EQ(2, p11[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p11[0]->elltwosquared());
    EXPECT_EQ(3, p11[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p11[1]->elltwosquared());
    EXPECT_EQ(5, p11[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p11[2]->elltwosquared());
    EXPECT_EQ(6, p11[3]->sum());
    EXPECT_EQ(4.5, p11[3]->elltwosquared());
    EXPECT_EQ(3, p11[4]->sum());
    EXPECT_EQ(4.5, p11[4]->elltwosquared());
    EXPECT_EQ(0, p11[6]->sum());
    EXPECT_EQ(4.5, p11[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.par_cbegin()}; it != e.par_cend(); ++it)
        ++i;
    EXPECT_EQ(e.par_size(), i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy assignment
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA>  p12(6, v.begin(), v.end());
    p12 = e;
    EXPECT_EQ(2ul, p12.dims());
    EXPECT_EQ(5ul, p12.par_size());
    EXPECT_EQ(6ul, p12.size());
    EXPECT_EQ(1ul, p12.count(0));
    EXPECT_EQ(1ul, p12.count(1));
    EXPECT_EQ(1ul, p12.count(2));
    EXPECT_EQ(1ul, p12.count(3));
    EXPECT_EQ(1ul, p12.count(4));
    EXPECT_EQ(0ul, p12.count(5));
    EXPECT_EQ(0ul, p12.count(6));
    i = 0;
    for (auto it {p12.cbegin()}; it != p12.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p12.refs(p12[0]));
    EXPECT_EQ(1ul, p12.refs(p12[1]));
    EXPECT_EQ(1ul, p12.refs(p12[2]));
    EXPECT_EQ(1ul, p12.refs(p12[3]));
    EXPECT_EQ(1ul, p12.refs(p12[4]));
    EXPECT_EQ(1ul, p12.refs(p12[6]));
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p12[0]->elltwosquared());
    EXPECT_EQ(3, p12[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p12[1]->elltwosquared());
    EXPECT_EQ(5, p12[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p12[2]->elltwosquared());
    EXPECT_EQ(6, p12[3]->sum());
    EXPECT_EQ(4.5, p12[3]->elltwosquared());
    EXPECT_EQ(3, p12[4]->sum());
    EXPECT_EQ(4.5, p12[4]->elltwosquared());
    EXPECT_EQ(0, p12[6]->sum());
    EXPECT_EQ(4.5, p12[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test move constructor
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p13 {std::move(e)};
    EXPECT_EQ(2ul, p13.dims());
    EXPECT_EQ(5ul, p13.par_size());
    EXPECT_EQ(6ul, p13.size());
    EXPECT_EQ(1ul, p13.count(0));
    EXPECT_EQ(1ul, p13.count(1));
    EXPECT_EQ(1ul, p13.count(2));
    EXPECT_EQ(1ul, p13.count(3));
    EXPECT_EQ(1ul, p13.count(4));
    EXPECT_EQ(0ul, p13.count(5));
    EXPECT_EQ(0ul, p13.count(6));
    i = 0;
    for (auto it {p13.cbegin()}; it != p13.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p13.refs(p13[0]));
    EXPECT_EQ(1ul, p13.refs(p13[1]));
    EXPECT_EQ(1ul, p13.refs(p13[2]));
    EXPECT_EQ(1ul, p13.refs(p13[3]));
    EXPECT_EQ(1ul, p13.refs(p13[4]));
    EXPECT_EQ(1ul, p13.refs(p13[6]));
    EXPECT_EQ(2, p13[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p13[0]->elltwosquared());
    EXPECT_EQ(3, p13[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p13[1]->elltwosquared());
    EXPECT_EQ(5, p13[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p13[2]->elltwosquared());
    EXPECT_EQ(6, p13[3]->sum());
    EXPECT_EQ(4.5, p13[3]->elltwosquared());
    EXPECT_EQ(3, p13[4]->sum());
    EXPECT_EQ(4.5, p13[4]->elltwosquared());
    EXPECT_EQ(0, p13[6]->sum());
    EXPECT_EQ(4.5, p13[6]->elltwosquared());
    // Test move assignment
    Wheelwright::Projection<int,Wheelwright::Features::DISCRETE | Wheelwright::Features::VANILLA> p14(6, v.begin(), v.end());
    p14 = std::move(p12);
    EXPECT_EQ(2ul, p14.dims());
    EXPECT_EQ(5ul, p14.par_size());
    EXPECT_EQ(6ul, p14.size());
    EXPECT_EQ(1ul, p14.count(0));
    EXPECT_EQ(1ul, p14.count(1));
    EXPECT_EQ(1ul, p14.count(2));
    EXPECT_EQ(1ul, p14.count(3));
    EXPECT_EQ(1ul, p14.count(4));
    EXPECT_EQ(0ul, p14.count(5));
    EXPECT_EQ(0ul, p14.count(6));
    i = 0;
    for (auto it {p14.cbegin()}; it != p14.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p14.refs(p14[0]));
    EXPECT_EQ(1ul, p14.refs(p14[1]));
    EXPECT_EQ(1ul, p14.refs(p14[2]));
    EXPECT_EQ(1ul, p14.refs(p14[3]));
    EXPECT_EQ(1ul, p14.refs(p14[4]));
    EXPECT_EQ(1ul, p14.refs(p14[6]));
    EXPECT_EQ(2, p14[0]->sum());
    EXPECT_DOUBLE_EQ(0.5, p14[0]->elltwosquared());
    EXPECT_EQ(3, p14[1]->sum());
    EXPECT_DOUBLE_EQ(0.5, p14[1]->elltwosquared());
    EXPECT_EQ(5, p14[2]->sum());
    EXPECT_DOUBLE_EQ(2.5, p14[2]->elltwosquared());
    EXPECT_EQ(6, p14[3]->sum());
    EXPECT_EQ(4.5, p14[3]->elltwosquared());
    EXPECT_EQ(3, p14[4]->sum());
    EXPECT_EQ(4.5, p14[4]->elltwosquared());
    EXPECT_EQ(0, p14[6]->sum());
    EXPECT_EQ(4.5, p14[6]->elltwosquared());
}

TEST(ProjectionTest, ProjectionWorksWithDiscrete) {
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    std::set<int> v {1, 2, 3};
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p(7,
            v.begin(), v.end());
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    EXPECT_EQ(nullptr, p[0]);
    // Test update()
    p.begin_update(0, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(2ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(6ul, p.refs(p[7]));
    EXPECT_NE(p[7]->hash(), p[0]->hash());
    std::size_t i {0};
    for (auto it {p.cbegin()}; it != p.cend(); ++it) {
        bool found_one {false};
        bool found_two {false};
        bool found_three {false};
        for (auto it2 {it->domain_begin()}; it2 != it->domain_end(); it2++) {
            switch (*it2) {
                case 1:
                    EXPECT_FALSE(found_one);
                    found_one = true;
                    break;
                case 2:
                    EXPECT_FALSE(found_two);
                    found_two = true;
                    break;
                case 3:
                    EXPECT_FALSE(found_three);
                    found_three = true;
                    break;
                default:
                    EXPECT_FALSE(true);
                    break;
            }
        }
        EXPECT_TRUE(found_one);
        EXPECT_TRUE(found_two);
        EXPECT_TRUE(found_three);
        ++i;
    }
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(0ul, p[0]->count(1));
    EXPECT_EQ(1ul, p[0]->count(2));
    EXPECT_EQ(0ul, p[0]->count(3));
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0ul, p[7]->count(1));
    EXPECT_EQ(0ul, p[7]->count(2));
    EXPECT_EQ(0ul, p[7]->count(3));
    // Test dims(d)
    p.dims(2);
    EXPECT_EQ(2ul, p.dims());
    // Test drop_catchall()
    p.drop_catchall();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(1ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(nullptr, p[7]);
    // Test clear()
    p.clear();
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    EXPECT_EQ(nullptr, p[0]);
    // Re-populate the Projection in preparation of testing update().
    p.begin_update(0, 1);
    p.update(0, 1);
    p.update(1, 1);
    p.update(2, 1);
    p.update(3, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(1ul, p[0]->count(1));
    EXPECT_EQ(0ul, p[0]->count(2));
    EXPECT_EQ(0ul, p[0]->count(3));
    EXPECT_EQ(2, p[3]->sum());
    EXPECT_EQ(0ul, p[3]->count(1));
    EXPECT_EQ(1ul, p[3]->count(2));
    EXPECT_EQ(0ul, p[3]->count(3));
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0ul, p[7]->count(1));
    EXPECT_EQ(0ul, p[7]->count(2));
    EXPECT_EQ(0ul, p[7]->count(3));
    // update() of parameter that is the only parameter mapped to its
    // vector
    p.begin_update(1, 1);
    p.update(3, 1);
    p.end_update();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(1ul, p[0]->count(1));
    EXPECT_EQ(0ul, p[0]->count(2));
    EXPECT_EQ(0ul, p[0]->count(3));
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_EQ(1ul, p[3]->count(1));
    EXPECT_EQ(1ul, p[3]->count(2));
    EXPECT_EQ(0ul, p[3]->count(3));
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0ul, p[7]->count(1));
    EXPECT_EQ(0ul, p[7]->count(2));
    EXPECT_EQ(0ul, p[7]->count(3));
    // update() with new parameter causing a new vector to be created.
    p.begin_update(2, 2);
    p.update(4, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(4ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(1ul, p[0]->count(1));
    EXPECT_EQ(0ul, p[0]->count(2));
    EXPECT_EQ(0ul, p[0]->count(3));
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_EQ(1ul, p[3]->count(1));
    EXPECT_EQ(1ul, p[3]->count(2));
    EXPECT_EQ(0ul, p[3]->count(3));
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_EQ(0ul, p[4]->count(1));
    EXPECT_EQ(2ul, p[4]->count(2));
    EXPECT_EQ(0ul, p[4]->count(3));
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0ul, p[7]->count(1));
    EXPECT_EQ(0ul, p[7]->count(2));
    EXPECT_EQ(0ul, p[7]->count(3));
    // Test update() with parameter which before the update is mapped to
    // a vector to which other parameters are mapped, and after the update is
    // mapped to a new vector.
    p.begin_update(3, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[1], p[3]);
    EXPECT_NE(p[1], p[4]);
    EXPECT_NE(p[1], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[1]->hash(), p[3]->hash());
    EXPECT_NE(p[1]->hash(), p[4]->hash());
    EXPECT_NE(p[1]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(3, p[0]->sum());
    EXPECT_EQ(1ul, p[0]->count(1));
    EXPECT_EQ(1ul, p[0]->count(2));
    EXPECT_EQ(0ul, p[0]->count(3));
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1ul, p[1]->count(1));
    EXPECT_EQ(0ul, p[1]->count(2));
    EXPECT_EQ(0ul, p[1]->count(3));
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_EQ(0ul, p[4]->count(1));
    EXPECT_EQ(2ul, p[4]->count(2));
    EXPECT_EQ(0ul, p[4]->count(3));
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0ul, p[7]->count(1));
    EXPECT_EQ(0ul, p[7]->count(2));
    EXPECT_EQ(0ul, p[7]->count(3));
    // Test merge() -- 1
    // Merge two Projections such that:
    //  * one has the catch-all vector and the other doesn't;
    //  * there is a parameter that is mapped to the catch-all vector in one of
    //    them, and to a non-catch-all vector in the other.
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p1(2,
            v.begin(), v.end());
    p1.begin_update(0, 1);
    p1.update(0, 1);
    p1.end_update();
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p2(2,
            v.begin(), v.end());
    p2.begin_update(1, 1);
    p2.update(0, 1);
    p2.update(1, 1);
    p2.end_update();
    p1.merge(p2);
    EXPECT_EQ(2ul, p1.dims());
    EXPECT_EQ(2ul, p1.par_size());
    EXPECT_EQ(2ul, p1.size());
    EXPECT_EQ(1ul, p1.count(0));
    EXPECT_EQ(1ul, p1.count(1));
    EXPECT_EQ(0ul, p1.count(2));
    i = 0;
    for (auto it {p1.cbegin()}; it != p1.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p1.refs(p1[0]));
    EXPECT_EQ(1ul, p1.refs(p1[1]));
    EXPECT_NE(p1[0], p1[1]);
    EXPECT_EQ(nullptr, p1[2]);
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(2, p1[0]->sum());
    EXPECT_EQ(1, p1[1]->sum());
    EXPECT_EQ(0ul, p2.dims());
    EXPECT_EQ(0ul, p2.par_size());
    EXPECT_EQ(0ul, p2.size());
    EXPECT_EQ(0ul, p2.count(0));
    EXPECT_EQ(0ul, p2.count(1));
    EXPECT_EQ(0ul, p2.count(2));
    i = 0;
    for (auto it {p2.cbegin()}; it != p2.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test merge() -- 2
    // Symmetric of the previous test.
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p3(2,
            v.begin(), v.end());
    p3.begin_update(1, 1);
    p3.update(0, 1);
    p3.update(1, 1);
    p3.end_update();
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p4(2,
            v.begin(), v.end());
    p4.begin_update(0, 1);
    p4.update(0, 1);
    p4.end_update();
    p3.merge(p4);
    EXPECT_EQ(2ul, p3.dims());
    EXPECT_EQ(2ul, p3.par_size());
    EXPECT_EQ(2ul, p3.size());
    EXPECT_EQ(1ul, p3.count(0));
    EXPECT_EQ(1ul, p3.count(1));
    EXPECT_EQ(0ul, p3.count(2));
    i = 0;
    for (auto it {p3.cbegin()}; it != p3.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p3.refs(p3[0]));
    EXPECT_EQ(1ul, p3.refs(p3[1]));
    EXPECT_NE(p3[0], p3[1]);
    EXPECT_EQ(nullptr, p3[2]);
    EXPECT_NE(p3[0]->hash(), p3[1]->hash());
    EXPECT_EQ(2, p3[0]->sum());
    EXPECT_EQ(1, p3[1]->sum());
    EXPECT_EQ(0ul, p4.dims());
    EXPECT_EQ(0ul, p4.par_size());
    EXPECT_EQ(0ul, p4.size());
    EXPECT_EQ(0ul, p4.count(0));
    EXPECT_EQ(0ul, p4.count(1));
    EXPECT_EQ(0ul, p4.count(2));
    i = 0;
    for (auto it {p4.cbegin()}; it != p4.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test merge() -- 3
    // Merge two Projections both containing the catch-all vector but such that
    // the merged one does not contain the catch-all vector.
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p5(2,
            v.begin(), v.end());
    p5.begin_update(0, 1);
    p5.update(0, 1);
    p5.end_update();
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p6(2,
            v.begin(), v.end());
    p6.begin_update(1, 1);
    p6.update(1, 1);
    p6.end_update();
    p5.merge(p6);
    EXPECT_EQ(2ul, p5.dims());
    EXPECT_EQ(2ul, p5.par_size());
    EXPECT_EQ(2ul, p5.size());
    EXPECT_EQ(1ul, p5.count(0));
    EXPECT_EQ(1ul, p5.count(1));
    EXPECT_EQ(0ul, p5.count(2));
    i = 0;
    for (auto it {p5.cbegin()}; it != p5.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p5.refs(p5[0]));
    EXPECT_EQ(1ul, p5.refs(p5[1]));
    EXPECT_NE(p5[0], p5[1]);
    EXPECT_EQ(nullptr, p5[2]);
    EXPECT_NE(p5[0]->hash(), p5[1]->hash());
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(1, p5[1]->sum());
    EXPECT_EQ(0ul, p6.dims());
    EXPECT_EQ(0ul, p6.par_size());
    EXPECT_EQ(0ul, p6.size());
    EXPECT_EQ(0ul, p6.count(0));
    EXPECT_EQ(0ul, p6.count(1));
    EXPECT_EQ(0ul, p6.count(2));
    i = 0;
    for (auto it {p6.cbegin()}; it != p6.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test merge() -- 4
    // Merge two Projections both containing the catch-all vector and such that
    // the merged one also contains the catch-all vector.
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p7(2,
            v.begin(), v.end());
    p7.begin_update(0, 1);
    p7.update(0, 1);
    p7.end_update();
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p8(2,
             v.begin(), v.end());
    p8.begin_update(1, 1);
    p8.update(0, 1);
    p8.end_update();
    p7.merge(p8);
    EXPECT_EQ(2ul, p7.dims());
    EXPECT_EQ(1ul, p7.par_size());
    EXPECT_EQ(2ul, p7.size());
    EXPECT_EQ(1ul, p7.count(0));
    EXPECT_EQ(0ul, p7.count(1));
    EXPECT_EQ(0ul, p7.count(2));
    i = 0;
    for (auto it {p7.cbegin()}; it != p7.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p7.refs(p7[0]));
    EXPECT_EQ(1ul, p7.refs(p7[2]));
    EXPECT_NE(p7[0], p7[2]);
    EXPECT_EQ(2, p7[0]->sum());
    EXPECT_EQ(0, p7[2]->sum());
    EXPECT_EQ(0ul, p8.dims());
    EXPECT_EQ(0ul, p8.par_size());
    EXPECT_EQ(0ul, p8.size());
    EXPECT_EQ(0ul, p8.count(0));
    EXPECT_EQ(0ul, p8.count(1));
    EXPECT_EQ(0ul, p8.count(2));
    i = 0;
    for (auto it {p8.cbegin()}; it != p8.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test merge() -- 5
    // Merge two Projections such that a pair of parameters is mapped to the
    // same vector in one Projection and to two different vectors in the other.
    // Vice versa for another pair of parameters.
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p9(6,
            v.begin(), v.end());
    p9.begin_update(0, 1);
    p9.update(0, 1);
    p9.update(1, 1);
    p9.update(2, 2);
    p9.update(3, 3);
    p9.update(4, 3);
    p9.end_update();
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p10(6,
            v.begin(), v.end());
    p10.begin_update(1, 1);
    p10.update(0, 1);
    p10.update(1, 2);
    p10.update(2, 3);
    p10.update(3, 3);
    p10.end_update();
    p9.merge(p10);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0ul, p10.dims());
    EXPECT_EQ(0ul, p10.par_size());
    EXPECT_EQ(0ul, p10.size());
    EXPECT_EQ(0ul, p10.count(0));
    EXPECT_EQ(0ul, p10.count(1));
    EXPECT_EQ(0ul, p10.count(2));
    i = 0;
    for (auto it {p10.cbegin()}; it != p10.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test merge()  -- 6
    // Test the special case that the oth Projection is empty.
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> e(6,
            v.begin(), v.end());
    p9.merge(e);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0ul, e.dims());
    EXPECT_EQ(0ul, e.par_size());
    EXPECT_EQ(0ul, e.size());
    EXPECT_EQ(0ul, e.count(0));
    EXPECT_EQ(0ul, e.count(1));
    EXPECT_EQ(0ul, e.count(2));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test merge() -- 7
    // Test the special case that this Projection is empty.
    e.merge(p9);
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(0, e[6]->sum());
    // Test copy constructor;
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p11 {e};
    EXPECT_EQ(2ul, p11.dims());
    EXPECT_EQ(5ul, p11.par_size());
    EXPECT_EQ(6ul, p11.size());
    EXPECT_EQ(1ul, p11.count(0));
    EXPECT_EQ(1ul, p11.count(1));
    EXPECT_EQ(1ul, p11.count(2));
    EXPECT_EQ(1ul, p11.count(3));
    EXPECT_EQ(1ul, p11.count(4));
    EXPECT_EQ(0ul, p11.count(5));
    EXPECT_EQ(0ul, p11.count(6));
    i = 0;
    for (auto it {p11.cbegin()}; it != p11.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p11.refs(p11[0]));
    EXPECT_EQ(1ul, p11.refs(p11[1]));
    EXPECT_EQ(1ul, p11.refs(p11[2]));
    EXPECT_EQ(1ul, p11.refs(p11[3]));
    EXPECT_EQ(1ul, p11.refs(p11[4]));
    EXPECT_EQ(1ul, p11.refs(p11[6]));
    EXPECT_EQ(2, p11[0]->sum());
    EXPECT_EQ(3, p11[1]->sum());
    EXPECT_EQ(5, p11[2]->sum());
    EXPECT_EQ(6, p11[3]->sum());
    EXPECT_EQ(3, p11[4]->sum());
    EXPECT_EQ(0, p11[6]->sum());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(0, e[6]->sum());
    // Test copy assignment
    Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p12(6,
            v.begin(), v.end());
    p12 = e;
    EXPECT_EQ(2ul, p12.dims());
    EXPECT_EQ(5ul, p12.par_size());
    EXPECT_EQ(6ul, p12.size());
    EXPECT_EQ(1ul, p12.count(0));
    EXPECT_EQ(1ul, p12.count(1));
    EXPECT_EQ(1ul, p12.count(2));
    EXPECT_EQ(1ul, p12.count(3));
    EXPECT_EQ(1ul, p12.count(4));
    EXPECT_EQ(0ul, p12.count(5));
    EXPECT_EQ(0ul, p12.count(6));
    i = 0;
    for (auto it {p12.cbegin()}; it != p12.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p12.refs(p12[0]));
    EXPECT_EQ(1ul, p12.refs(p12[1]));
    EXPECT_EQ(1ul, p12.refs(p12[2]));
    EXPECT_EQ(1ul, p12.refs(p12[3]));
    EXPECT_EQ(1ul, p12.refs(p12[4]));
    EXPECT_EQ(1ul, p12.refs(p12[6]));
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_EQ(3, p12[1]->sum());
    EXPECT_EQ(5, p12[2]->sum());
    EXPECT_EQ(6, p12[3]->sum());
    EXPECT_EQ(3, p12[4]->sum());
    EXPECT_EQ(0, p12[6]->sum());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(0, e[6]->sum());
    // Test move constructor
     Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p13
     {std::move(e)};
    EXPECT_EQ(2ul, p13.dims());
    EXPECT_EQ(5ul, p13.par_size());
    EXPECT_EQ(6ul, p13.size());
    EXPECT_EQ(1ul, p13.count(0));
    EXPECT_EQ(1ul, p13.count(1));
    EXPECT_EQ(1ul, p13.count(2));
    EXPECT_EQ(1ul, p13.count(3));
    EXPECT_EQ(1ul, p13.count(4));
    EXPECT_EQ(0ul, p13.count(5));
    EXPECT_EQ(0ul, p13.count(6));
    i = 0;
    for (auto it {p13.cbegin()}; it != p13.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p13.refs(p13[0]));
    EXPECT_EQ(1ul, p13.refs(p13[1]));
    EXPECT_EQ(1ul, p13.refs(p13[2]));
    EXPECT_EQ(1ul, p13.refs(p13[3]));
    EXPECT_EQ(1ul, p13.refs(p13[4]));
    EXPECT_EQ(1ul, p13.refs(p13[6]));
    EXPECT_EQ(2, p13[0]->sum());
    EXPECT_EQ(3, p13[1]->sum());
    EXPECT_EQ(5, p13[2]->sum());
    EXPECT_EQ(6, p13[3]->sum());
    EXPECT_EQ(3, p13[4]->sum());
    EXPECT_EQ(0, p13[6]->sum());
    // Test move assignment
     Wheelwright::Projection<int, Wheelwright::Features::DISCRETE, int> p14(6,
             v.begin(), v.end());
    p14 = std::move(p12);
    EXPECT_EQ(2ul, p14.dims());
    EXPECT_EQ(5ul, p14.par_size());
    EXPECT_EQ(6ul, p14.size());
    EXPECT_EQ(1ul, p14.count(0));
    EXPECT_EQ(1ul, p14.count(1));
    EXPECT_EQ(1ul, p14.count(2));
    EXPECT_EQ(1ul, p14.count(3));
    EXPECT_EQ(1ul, p14.count(4));
    EXPECT_EQ(0ul, p14.count(5));
    EXPECT_EQ(0ul, p14.count(6));
    i = 0;
    for (auto it {p14.cbegin()}; it != p14.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p14.refs(p14[0]));
    EXPECT_EQ(1ul, p14.refs(p14[1]));
    EXPECT_EQ(1ul, p14.refs(p14[2]));
    EXPECT_EQ(1ul, p14.refs(p14[3]));
    EXPECT_EQ(1ul, p14.refs(p14[4]));
    EXPECT_EQ(1ul, p14.refs(p14[6]));
    EXPECT_EQ(2, p14[0]->sum());
    EXPECT_EQ(3, p14[1]->sum());
    EXPECT_EQ(5, p14[2]->sum());
    EXPECT_EQ(6, p14[3]->sum());
    EXPECT_EQ(3, p14[4]->sum());
    EXPECT_EQ(0, p14[6]->sum());
}

TEST(ProjectionTest, ProjectionWorksWithOneDraw) {
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    Wheelwright::Projection<int, Wheelwright::Features::ONEDRAW, int> p(7);
    //set the seed for reproducibility.
    p.seed(1);
    // We can use the onedrawsum() function of the vectors after each call to
    // p.end_update() to test that the Rademacher averages have been set
    // correctly.
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    EXPECT_EQ(nullptr, p[0]);
    std::size_t i {0};
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test populating the Projection (simple case to then test drop_catchall
    // and clear. A more complex test is done later)
    p.begin_update(0, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(2ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(6ul, p.refs(p[7]));
    EXPECT_NE(p[7]->hash(), p[0]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    // Test dims(d)
    p.dims(2);
    EXPECT_EQ(2ul, p.dims());
    // Test drop_catchall()
    p.drop_catchall();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(1ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(nullptr, p[7]);
    // Test clear()
    p.clear();
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    EXPECT_EQ(nullptr, p[0]);
    // Test starting to populate the Projection using the
    // begin_update()-update()-end_update() sequence.
    p.begin_update(0, 1);
    p.update(0, 1);
    p.update(1, 1);
    p.update(2, 1);
    p.update(3, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(2, p[3]->sum());
    EXPECT_EQ(2, p[3]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    // update() of parameter that is the only parameter mapped to its vector
    p.begin_update(1, 1);
    p.update(3, 1);
    p.end_update();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    // update() with new parameter causing a new vector to be created.
    p.begin_update(2, 2);
    p.update(4, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(4ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    // Test update() with parameter which before the update is mapped to a
    // vector to which other parameters are mapped, and after the update is
    // mapped to a new vector.
    p.begin_update(4, 2);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(6ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[1], p[3]);
    EXPECT_NE(p[1], p[4]);
    EXPECT_NE(p[1], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[1]->hash(), p[3]->hash());
    EXPECT_NE(p[1]->hash(), p[4]->hash());
    EXPECT_NE(p[1]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    // Test centralized update dropping the cav vector when needed.
    p.begin_update(6, 1);
    p.update(5, 5);
    p.update(6, 5);
    p.end_update();
    EXPECT_EQ(7ul, p.dims());
    EXPECT_EQ(7ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(1ul, p.count(5));
    EXPECT_EQ(1ul, p.count(6));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[5]));
    EXPECT_EQ(2ul, p.refs(p[6]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[5]);
    EXPECT_NE(p[0], p[6]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[2], p[3]);
    EXPECT_NE(p[2], p[4]);
    EXPECT_NE(p[2], p[5]);
    EXPECT_NE(p[2], p[6]);
    EXPECT_NE(p[2], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[3], p[5]);
    EXPECT_NE(p[3], p[6]);
    EXPECT_NE(p[3], p[7]);
    EXPECT_NE(p[4], p[5]);
    EXPECT_NE(p[4], p[6]);
    EXPECT_NE(p[4], p[7]);
    EXPECT_EQ(p[5], p[6]);
    EXPECT_NE(p[5], p[7]);
    EXPECT_EQ(nullptr, p[7]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[5]->hash());
    EXPECT_NE(p[0]->hash(), p[6]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[2]->hash(), p[3]->hash());
    EXPECT_NE(p[2]->hash(), p[4]->hash());
    EXPECT_NE(p[2]->hash(), p[5]->hash());
    EXPECT_NE(p[2]->hash(), p[6]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[3]->hash(), p[5]->hash());
    EXPECT_NE(p[3]->hash(), p[6]->hash());
    EXPECT_NE(p[4]->hash(), p[5]->hash());
    EXPECT_NE(p[4]->hash(), p[6]->hash());
    EXPECT_EQ(p[5]->hash(), p[6]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_EQ(5, p[5]->sum());
    EXPECT_DOUBLE_EQ(5, p[5]->onedrawsum());
    EXPECT_EQ(5, p[6]->sum());
    EXPECT_DOUBLE_EQ(5, p[6]->onedrawsum());
    // Test merge()  -- 1
    // Merge two Projections such that:
    //  * one has the catch-all vector and the other doesn't;
    //  * there is a parameter that is mapped to the catch-all vector in one of
    //    them, and to a non-catch-all vector in the other.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p1(2);
    p1.seed(1);
    p1.begin_update(0, 1);
    p1.update(0, 1);
    p1.end_update();
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(1, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p2(2);
    p2.seed(2);
    p2.begin_update(1, 1);
    p2.update(0, 1);
    p2.update(1, 1);
    p2.end_update();
    EXPECT_EQ(1, p2[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[0]->onedrawsum());
    EXPECT_EQ(1, p2[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[1]->onedrawsum());
    p1.merge(p2);
    EXPECT_EQ(2ul, p1.dims());
    EXPECT_EQ(2ul, p1.par_size());
    EXPECT_EQ(2ul, p1.size());
    EXPECT_EQ(1ul, p1.count(0));
    EXPECT_EQ(1ul, p1.count(1));
    EXPECT_EQ(0ul, p1.count(2));
    i = 0;
    for (auto it {p1.cbegin()}; it != p1.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p1.refs(p1[0]));
    EXPECT_EQ(1ul, p1.refs(p1[1]));
    EXPECT_NE(p1[0], p1[1]);
    EXPECT_EQ(nullptr, p1[2]);
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(2, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p1[0]->onedrawsum());
    EXPECT_EQ(1, p1[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[1]->onedrawsum());
    // Test merge() -- 2
    // Symmetric of the previous test.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p3(2);
    p3.seed(2);
    p3.begin_update(1, 1);
    p3.update(0, 1);
    p3.update(1, 1);
    p3.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p4(2);
    p4.seed(1);
    p4.begin_update(0, 1);
    p4.update(0, 1);
    p4.end_update();
    p3.merge(p4);
    EXPECT_EQ(2ul, p3.dims());
    EXPECT_EQ(2ul, p3.par_size());
    EXPECT_EQ(2ul, p3.size());
    EXPECT_EQ(1ul, p3.count(0));
    EXPECT_EQ(1ul, p3.count(1));
    EXPECT_EQ(0ul, p3.count(2));
    i = 0;
    for (auto it {p3.cbegin()}; it != p3.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p3.refs(p3[0]));
    EXPECT_EQ(1ul, p3.refs(p3[1]));
    EXPECT_NE(p3[0], p3[1]);
    EXPECT_EQ(nullptr, p3[2]);
    EXPECT_NE(p3[0]->hash(), p3[1]->hash());
    EXPECT_EQ(2, p3[0]->sum());
    EXPECT_EQ(-2, p3[0]->onedrawsum());
    EXPECT_EQ(1, p3[1]->sum());
    EXPECT_EQ(-1, p3[1]->onedrawsum());
    // Test merge() -- 3
    // Merge two Projections both containing the catch-all vector but such that
    // the merged one does not contain the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p5(2);
    p5.seed(1);
    p5.begin_update(0, 1);
    p5.update(0, 1);
    p5.end_update();
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p6(2);
    p6.seed(2);
    p6.begin_update(1, 1);
    p6.update(1, 1);
    p6.end_update();
    EXPECT_EQ(1, p6[1]->sum());
    EXPECT_EQ(-1, p6[1]->onedrawsum());
    p5.merge(p6);
    EXPECT_EQ(2ul, p5.dims());
    EXPECT_EQ(2ul, p5.par_size());
    EXPECT_EQ(2ul, p5.size());
    EXPECT_EQ(1ul, p5.count(0));
    EXPECT_EQ(1ul, p5.count(1));
    EXPECT_EQ(0ul, p5.count(2));
    i = 0;
    for (auto it {p5.cbegin()}; it != p5.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p5.refs(p5[0]));
    EXPECT_EQ(1ul, p5.refs(p5[1]));
    EXPECT_NE(p5[0], p5[1]);
    EXPECT_EQ(nullptr, p5[2]);
    EXPECT_NE(p5[0]->hash(), p5[1]->hash());
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    EXPECT_EQ(1, p5[1]->sum());
    EXPECT_EQ(-1, p5[1]->onedrawsum());
    // Test merge() -- 4
    // Merge two Projections both containing the catch-all vector and such that
    // the merged one also contains the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p7(2);
    p7.seed(1);
    p7.begin_update(0, 1);
    p7.update(0, 1);
    p7.end_update();
    EXPECT_EQ(1, p7[0]->sum());
    EXPECT_EQ(-1, p7[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p8(2);
    p8.seed(2);
    p8.begin_update(1, 1);
    p8.update(0, 1);
    p8.end_update();
    EXPECT_EQ(1, p8[0]->sum());
    EXPECT_EQ(-1, p8[0]->onedrawsum());
    p7.merge(p8);
    EXPECT_EQ(2ul, p7.dims());
    EXPECT_EQ(1ul, p7.par_size());
    EXPECT_EQ(2ul, p7.size());
    EXPECT_EQ(1ul, p7.count(0));
    EXPECT_EQ(0ul, p7.count(1));
    EXPECT_EQ(0ul, p7.count(2));
    i = 0;
    for (auto it {p7.cbegin()}; it != p7.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p7.refs(p7[0]));
    EXPECT_EQ(1ul, p7.refs(p7[2]));
    EXPECT_NE(p7[0], p7[2]);
    EXPECT_EQ(2, p7[0]->sum());
    EXPECT_EQ(-2, p7[0]->onedrawsum());
    EXPECT_EQ(0, p7[2]->sum());
    EXPECT_EQ(0, p7[2]->onedrawsum());
    // Test merge()  -- 5
    // Merge two Projections such that a pair of parameters is mapped to the
    // same vector in one Projection and to two different vectors in the other.
    // Vice versa for another pair of parameters.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p9(6);
    p9.seed(1);
    p9.begin_update(0, 1);
    p9.update(0, 1);
    p9.update(1, 1);
    p9.update(2, 2);
    p9.update(3, 3);
    p9.update(4, 3);
    p9.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p10(6);
    p10.seed(2);
    p10.begin_update(1, 1);
    p10.update(0, 1);
    p10.update(1, 2);
    p10.update(2, 3);
    p10.update(3, 3);
    p10.end_update();
    p9.merge(p10);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    // Test merge()  -- 6
    // Test the special case that the oth Projection is empty.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> e(6);
    p9.merge(e);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    // Test merge() -- 7
    // Test the special case that this Projection is empty.
    e.merge(p9);
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    // Test copy constructor;
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p11 {e};
    EXPECT_EQ(2ul, p11.dims());
    EXPECT_EQ(5ul, p11.par_size());
    EXPECT_EQ(6ul, p11.size());
    EXPECT_EQ(1ul, p11.count(0));
    EXPECT_EQ(1ul, p11.count(1));
    EXPECT_EQ(1ul, p11.count(2));
    EXPECT_EQ(1ul, p11.count(3));
    EXPECT_EQ(1ul, p11.count(4));
    EXPECT_EQ(0ul, p11.count(5));
    EXPECT_EQ(0ul, p11.count(6));
    i = 0;
    for (auto it {p11.cbegin()}; it != p11.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p11.refs(p11[0]));
    EXPECT_EQ(1ul, p11.refs(p11[1]));
    EXPECT_EQ(1ul, p11.refs(p11[2]));
    EXPECT_EQ(1ul, p11.refs(p11[3]));
    EXPECT_EQ(1ul, p11.refs(p11[4]));
    EXPECT_EQ(1ul, p11.refs(p11[6]));
    EXPECT_EQ(2, p11[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p11[0]->onedrawsum());
    EXPECT_EQ(3, p11[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p11[1]->onedrawsum());
    EXPECT_EQ(5, p11[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p11[2]->onedrawsum());
    EXPECT_EQ(6, p11[3]->sum());
    EXPECT_EQ(-6, p11[3]->onedrawsum());
    EXPECT_EQ(3, p11[4]->sum());
    EXPECT_EQ(-3, p11[4]->onedrawsum());
    EXPECT_EQ(0, p11[6]->sum());
    EXPECT_EQ(0, p11[6]->onedrawsum());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.par_cbegin()}; it != e.par_cend(); ++it)
        ++i;
    EXPECT_EQ(e.par_size(), i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    // Test copy assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW>  p12(6);
    p12 = e;
    EXPECT_EQ(2ul, p12.dims());
    EXPECT_EQ(5ul, p12.par_size());
    EXPECT_EQ(6ul, p12.size());
    EXPECT_EQ(1ul, p12.count(0));
    EXPECT_EQ(1ul, p12.count(1));
    EXPECT_EQ(1ul, p12.count(2));
    EXPECT_EQ(1ul, p12.count(3));
    EXPECT_EQ(1ul, p12.count(4));
    EXPECT_EQ(0ul, p12.count(5));
    EXPECT_EQ(0ul, p12.count(6));
    i = 0;
    for (auto it {p12.cbegin()}; it != p12.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p12.refs(p12[0]));
    EXPECT_EQ(1ul, p12.refs(p12[1]));
    EXPECT_EQ(1ul, p12.refs(p12[2]));
    EXPECT_EQ(1ul, p12.refs(p12[3]));
    EXPECT_EQ(1ul, p12.refs(p12[4]));
    EXPECT_EQ(1ul, p12.refs(p12[6]));
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p12[0]->onedrawsum());
    EXPECT_EQ(3, p12[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p12[1]->onedrawsum());
    EXPECT_EQ(5, p12[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p12[2]->onedrawsum());
    EXPECT_EQ(6, p12[3]->sum());
    EXPECT_EQ(-6, p12[3]->onedrawsum());
    EXPECT_EQ(3, p12[4]->sum());
    EXPECT_EQ(-3, p12[4]->onedrawsum());
    EXPECT_EQ(0, p12[6]->sum());
    EXPECT_EQ(0, p12[6]->onedrawsum());
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    EXPECT_EQ(2, e[0]->sum());
    // Test move constructor
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p13 {std::move(e)};
    EXPECT_EQ(2ul, p13.dims());
    EXPECT_EQ(5ul, p13.par_size());
    EXPECT_EQ(6ul, p13.size());
    EXPECT_EQ(1ul, p13.count(0));
    EXPECT_EQ(1ul, p13.count(1));
    EXPECT_EQ(1ul, p13.count(2));
    EXPECT_EQ(1ul, p13.count(3));
    EXPECT_EQ(1ul, p13.count(4));
    EXPECT_EQ(0ul, p13.count(5));
    EXPECT_EQ(0ul, p13.count(6));
    i = 0;
    for (auto it {p13.cbegin()}; it != p13.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p13.refs(p13[0]));
    EXPECT_EQ(1ul, p13.refs(p13[1]));
    EXPECT_EQ(1ul, p13.refs(p13[2]));
    EXPECT_EQ(1ul, p13.refs(p13[3]));
    EXPECT_EQ(1ul, p13.refs(p13[4]));
    EXPECT_EQ(1ul, p13.refs(p13[6]));
    EXPECT_EQ(2, p13[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p13[0]->onedrawsum());
    EXPECT_EQ(3, p13[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p13[1]->onedrawsum());
    EXPECT_EQ(5, p13[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p13[2]->onedrawsum());
    EXPECT_EQ(6, p13[3]->sum());
    EXPECT_EQ(-6, p13[3]->onedrawsum());
    EXPECT_EQ(3, p13[4]->sum());
    EXPECT_EQ(-3, p13[4]->onedrawsum());
    EXPECT_EQ(0, p13[6]->sum());
    EXPECT_EQ(0, p13[6]->onedrawsum());
    EXPECT_EQ(2, p13[0]->sum());
    // Test move assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW> p14(6);
    p14 = std::move(p12);
    EXPECT_EQ(2ul, p14.dims());
    EXPECT_EQ(5ul, p14.par_size());
    EXPECT_EQ(6ul, p14.size());
    EXPECT_EQ(1ul, p14.count(0));
    EXPECT_EQ(1ul, p14.count(1));
    EXPECT_EQ(1ul, p14.count(2));
    EXPECT_EQ(1ul, p14.count(3));
    EXPECT_EQ(1ul, p14.count(4));
    EXPECT_EQ(0ul, p14.count(5));
    EXPECT_EQ(0ul, p14.count(6));
    i = 0;
    for (auto it {p14.cbegin()}; it != p14.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p14.refs(p14[0]));
    EXPECT_EQ(1ul, p14.refs(p14[1]));
    EXPECT_EQ(1ul, p14.refs(p14[2]));
    EXPECT_EQ(1ul, p14.refs(p14[3]));
    EXPECT_EQ(1ul, p14.refs(p14[4]));
    EXPECT_EQ(1ul, p14.refs(p14[6]));
    EXPECT_EQ(2, p14[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p14[0]->onedrawsum());
    EXPECT_EQ(3, p14[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p14[1]->onedrawsum());
    EXPECT_EQ(5, p14[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p14[2]->onedrawsum());
    EXPECT_EQ(6, p14[3]->sum());
    EXPECT_EQ(-6, p14[3]->onedrawsum());
    EXPECT_EQ(3, p14[4]->sum());
    EXPECT_EQ(-3, p14[4]->onedrawsum());
    EXPECT_EQ(0, p14[6]->sum());
    EXPECT_EQ(0, p14[6]->onedrawsum());
}


TEST(ProjectionTest, ProjectionWorksWithOneDrawAndDiscrete) {
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    std::set<int> v {1, 2, 3, 4, 5};
    Wheelwright::Projection<int, Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE, int> p(7, v.begin(), v.end());
    //set the seed for reproducibility.
    p.seed(1);
    // We can use the onedrawsum() function of the vectors after each call to
    // p.end_update() to test that the Rademacher averages have been set
    // correctly.
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    EXPECT_EQ(nullptr, p[0]);
    std::size_t i {0};
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test populating the Projection (simple case to then test drop_catchall
    // and clear. A more complex test is done later)
    p.begin_update(0, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(2ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(6ul, p.refs(p[7]));
    EXPECT_NE(p[7]->hash(), p[0]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it) {
        bool found_one {false};
        bool found_two {false};
        bool found_three {false};
        bool found_four {false};
        bool found_five {false};
        for (auto it2 {it->domain_begin()}; it2 != it->domain_end(); it2++) {
            switch (*it2) {
                case 1:
                    EXPECT_FALSE(found_one);
                    found_one = true;
                    break;
                case 2:
                    EXPECT_FALSE(found_two);
                    found_two = true;
                    break;
                case 3:
                    EXPECT_FALSE(found_three);
                    found_three = true;
                    break;
                case 4:
                    EXPECT_FALSE(found_four);
                    found_four = true;
                    break;
                case 5:
                    EXPECT_FALSE(found_five);
                    found_five = true;
                    break;
                default:
                    EXPECT_FALSE(true);
                    break;
            }
        }
        EXPECT_TRUE(found_one);
        EXPECT_TRUE(found_two);
        EXPECT_TRUE(found_three);
        EXPECT_TRUE(found_four);
        EXPECT_TRUE(found_five);
        ++i;
    }
    EXPECT_EQ(2ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    // Test dims(d)
    p.dims(2);
    EXPECT_EQ(2ul, p.dims());
    // Test drop_catchall()
    p.drop_catchall();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(1ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(nullptr, p[7]);
    // Test clear()
    p.clear();
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    EXPECT_EQ(nullptr, p[0]);
    // Test starting to populate the Projection using the
    // begin_update()-update()-end_update() sequence.
    p.begin_update(0, 1);
    p.update(0, 1);
    p.update(1, 1);
    p.update(2, 1);
    p.update(3, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(2, p[3]->sum());
    EXPECT_EQ(2, p[3]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    // update() of parameter that is the only parameter mapped to its vector
    p.begin_update(1, 1);
    p.update(3, 1);
    p.end_update();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    // update() with new parameter causing a new vector to be created.
    p.begin_update(2, 2);
    p.update(4, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(4ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    // Test update() with parameter which before the update is mapped to a
    // vector to which other parameters are mapped, and after the update is
    // mapped to a new vector.
    p.begin_update(4, 2);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(6ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[1], p[3]);
    EXPECT_NE(p[1], p[4]);
    EXPECT_NE(p[1], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[1]->hash(), p[3]->hash());
    EXPECT_NE(p[1]->hash(), p[4]->hash());
    EXPECT_NE(p[1]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    // Test centralized update dropping the cav vector when needed.
    p.begin_update(6, 1);
    p.update(5, 5);
    p.update(6, 5);
    p.end_update();
    EXPECT_EQ(7ul, p.dims());
    EXPECT_EQ(7ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(1ul, p.count(5));
    EXPECT_EQ(1ul, p.count(6));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[5]));
    EXPECT_EQ(2ul, p.refs(p[6]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[5]);
    EXPECT_NE(p[0], p[6]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[2], p[3]);
    EXPECT_NE(p[2], p[4]);
    EXPECT_NE(p[2], p[5]);
    EXPECT_NE(p[2], p[6]);
    EXPECT_NE(p[2], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[3], p[5]);
    EXPECT_NE(p[3], p[6]);
    EXPECT_NE(p[3], p[7]);
    EXPECT_NE(p[4], p[5]);
    EXPECT_NE(p[4], p[6]);
    EXPECT_NE(p[4], p[7]);
    EXPECT_EQ(p[5], p[6]);
    EXPECT_NE(p[5], p[7]);
    EXPECT_EQ(nullptr, p[7]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[5]->hash());
    EXPECT_NE(p[0]->hash(), p[6]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[2]->hash(), p[3]->hash());
    EXPECT_NE(p[2]->hash(), p[4]->hash());
    EXPECT_NE(p[2]->hash(), p[5]->hash());
    EXPECT_NE(p[2]->hash(), p[6]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[3]->hash(), p[5]->hash());
    EXPECT_NE(p[3]->hash(), p[6]->hash());
    EXPECT_NE(p[4]->hash(), p[5]->hash());
    EXPECT_NE(p[4]->hash(), p[6]->hash());
    EXPECT_EQ(p[5]->hash(), p[6]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_EQ(5, p[5]->sum());
    EXPECT_DOUBLE_EQ(5, p[5]->onedrawsum());
    EXPECT_EQ(5, p[6]->sum());
    EXPECT_DOUBLE_EQ(5, p[6]->onedrawsum());
    // Test merge()  -- 1
    // Merge two Projections such that:
    //  * one has the catch-all vector and the other doesn't;
    //  * there is a parameter that is mapped to the catch-all vector in one of
    //    them, and to a non-catch-all vector in the other.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p1(2, v.begin(), v.end());
    p1.seed(1);
    p1.begin_update(0, 1);
    p1.update(0, 1);
    p1.end_update();
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(1, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p2(2, v.begin(), v.end());
    p2.seed(2);
    p2.begin_update(1, 1);
    p2.update(0, 1);
    p2.update(1, 1);
    p2.end_update();
    EXPECT_EQ(1, p2[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[0]->onedrawsum());
    EXPECT_EQ(1, p2[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[1]->onedrawsum());
    p1.merge(p2);
    EXPECT_EQ(2ul, p1.dims());
    EXPECT_EQ(2ul, p1.par_size());
    EXPECT_EQ(2ul, p1.size());
    EXPECT_EQ(1ul, p1.count(0));
    EXPECT_EQ(1ul, p1.count(1));
    EXPECT_EQ(0ul, p1.count(2));
    i = 0;
    for (auto it {p1.cbegin()}; it != p1.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p1.refs(p1[0]));
    EXPECT_EQ(1ul, p1.refs(p1[1]));
    EXPECT_NE(p1[0], p1[1]);
    EXPECT_EQ(nullptr, p1[2]);
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(2, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p1[0]->onedrawsum());
    EXPECT_EQ(1, p1[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[1]->onedrawsum());
    // Test merge() -- 2
    // Symmetric of the previous test.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p3(2, v.begin(), v.end());
    p3.seed(2);
    p3.begin_update(1, 1);
    p3.update(0, 1);
    p3.update(1, 1);
    p3.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p4(2, v.begin(), v.end());
    p4.seed(1);
    p4.begin_update(0, 1);
    p4.update(0, 1);
    p4.end_update();
    p3.merge(p4);
    EXPECT_EQ(2ul, p3.dims());
    EXPECT_EQ(2ul, p3.par_size());
    EXPECT_EQ(2ul, p3.size());
    EXPECT_EQ(1ul, p3.count(0));
    EXPECT_EQ(1ul, p3.count(1));
    EXPECT_EQ(0ul, p3.count(2));
    i = 0;
    for (auto it {p3.cbegin()}; it != p3.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p3.refs(p3[0]));
    EXPECT_EQ(1ul, p3.refs(p3[1]));
    EXPECT_NE(p3[0], p3[1]);
    EXPECT_EQ(nullptr, p3[2]);
    EXPECT_NE(p3[0]->hash(), p3[1]->hash());
    EXPECT_EQ(2, p3[0]->sum());
    EXPECT_EQ(-2, p3[0]->onedrawsum());
    EXPECT_EQ(1, p3[1]->sum());
    EXPECT_EQ(-1, p3[1]->onedrawsum());
    // Test merge() -- 3
    // Merge two Projections both containing the catch-all vector but such that
    // the merged one does not contain the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p5(2, v.begin(), v.end());
    p5.seed(1);
    p5.begin_update(0, 1);
    p5.update(0, 1);
    p5.end_update();
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p6(2, v.begin(), v.end());
    p6.seed(2);
    p6.begin_update(1, 1);
    p6.update(1, 1);
    p6.end_update();
    EXPECT_EQ(1, p6[1]->sum());
    EXPECT_EQ(-1, p6[1]->onedrawsum());
    p5.merge(p6);
    EXPECT_EQ(2ul, p5.dims());
    EXPECT_EQ(2ul, p5.par_size());
    EXPECT_EQ(2ul, p5.size());
    EXPECT_EQ(1ul, p5.count(0));
    EXPECT_EQ(1ul, p5.count(1));
    EXPECT_EQ(0ul, p5.count(2));
    i = 0;
    for (auto it {p5.cbegin()}; it != p5.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p5.refs(p5[0]));
    EXPECT_EQ(1ul, p5.refs(p5[1]));
    EXPECT_NE(p5[0], p5[1]);
    EXPECT_EQ(nullptr, p5[2]);
    EXPECT_NE(p5[0]->hash(), p5[1]->hash());
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    EXPECT_EQ(1, p5[1]->sum());
    EXPECT_EQ(-1, p5[1]->onedrawsum());
    // Test merge() -- 4
    // Merge two Projections both containing the catch-all vector and such that
    // the merged one also contains the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p7(2, v.begin(), v.end());
    p7.seed(1);
    p7.begin_update(0, 1);
    p7.update(0, 1);
    p7.end_update();
    EXPECT_EQ(1, p7[0]->sum());
    EXPECT_EQ(-1, p7[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p8(2, v.begin(), v.end());
    p8.seed(2);
    p8.begin_update(1, 1);
    p8.update(0, 1);
    p8.end_update();
    EXPECT_EQ(1, p8[0]->sum());
    EXPECT_EQ(-1, p8[0]->onedrawsum());
    p7.merge(p8);
    EXPECT_EQ(2ul, p7.dims());
    EXPECT_EQ(1ul, p7.par_size());
    EXPECT_EQ(2ul, p7.size());
    EXPECT_EQ(1ul, p7.count(0));
    EXPECT_EQ(0ul, p7.count(1));
    EXPECT_EQ(0ul, p7.count(2));
    i = 0;
    for (auto it {p7.cbegin()}; it != p7.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p7.refs(p7[0]));
    EXPECT_EQ(1ul, p7.refs(p7[2]));
    EXPECT_NE(p7[0], p7[2]);
    EXPECT_EQ(2, p7[0]->sum());
    EXPECT_EQ(-2, p7[0]->onedrawsum());
    EXPECT_EQ(0, p7[2]->sum());
    EXPECT_EQ(0, p7[2]->onedrawsum());
    // Test merge()  -- 5
    // Merge two Projections such that a pair of parameters is mapped to the
    // same vector in one Projection and to two different vectors in the other.
    // Vice versa for another pair of parameters.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p9(6, v.begin(), v.end());
    p9.seed(1);
    p9.begin_update(0, 1);
    p9.update(0, 1);
    p9.update(1, 1);
    p9.update(2, 2);
    p9.update(3, 3);
    p9.update(4, 3);
    p9.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p10(6, v.begin(), v.end());
    p10.seed(2);
    p10.begin_update(1, 1);
    p10.update(0, 1);
    p10.update(1, 2);
    p10.update(2, 3);
    p10.update(3, 3);
    p10.end_update();
    p9.merge(p10);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    // Test merge()  -- 6
    // Test the special case that the oth Projection is empty.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> e(6, v.begin(), v.end());
    p9.merge(e);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    // Test merge() -- 7
    // Test the special case that this Projection is empty.
    e.merge(p9);
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    // Test copy constructor;
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p11 {e};
    EXPECT_EQ(2ul, p11.dims());
    EXPECT_EQ(5ul, p11.par_size());
    EXPECT_EQ(6ul, p11.size());
    EXPECT_EQ(1ul, p11.count(0));
    EXPECT_EQ(1ul, p11.count(1));
    EXPECT_EQ(1ul, p11.count(2));
    EXPECT_EQ(1ul, p11.count(3));
    EXPECT_EQ(1ul, p11.count(4));
    EXPECT_EQ(0ul, p11.count(5));
    EXPECT_EQ(0ul, p11.count(6));
    i = 0;
    for (auto it {p11.cbegin()}; it != p11.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p11.refs(p11[0]));
    EXPECT_EQ(1ul, p11.refs(p11[1]));
    EXPECT_EQ(1ul, p11.refs(p11[2]));
    EXPECT_EQ(1ul, p11.refs(p11[3]));
    EXPECT_EQ(1ul, p11.refs(p11[4]));
    EXPECT_EQ(1ul, p11.refs(p11[6]));
    EXPECT_EQ(2, p11[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p11[0]->onedrawsum());
    EXPECT_EQ(3, p11[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p11[1]->onedrawsum());
    EXPECT_EQ(5, p11[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p11[2]->onedrawsum());
    EXPECT_EQ(6, p11[3]->sum());
    EXPECT_EQ(-6, p11[3]->onedrawsum());
    EXPECT_EQ(3, p11[4]->sum());
    EXPECT_EQ(-3, p11[4]->onedrawsum());
    EXPECT_EQ(0, p11[6]->sum());
    EXPECT_EQ(0, p11[6]->onedrawsum());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.par_cbegin()}; it != e.par_cend(); ++it)
        ++i;
    EXPECT_EQ(e.par_size(), i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    // Test copy assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE>  p12(6, v.begin(), v.end());
    p12 = e;
    EXPECT_EQ(2ul, p12.dims());
    EXPECT_EQ(5ul, p12.par_size());
    EXPECT_EQ(6ul, p12.size());
    EXPECT_EQ(1ul, p12.count(0));
    EXPECT_EQ(1ul, p12.count(1));
    EXPECT_EQ(1ul, p12.count(2));
    EXPECT_EQ(1ul, p12.count(3));
    EXPECT_EQ(1ul, p12.count(4));
    EXPECT_EQ(0ul, p12.count(5));
    EXPECT_EQ(0ul, p12.count(6));
    i = 0;
    for (auto it {p12.cbegin()}; it != p12.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p12.refs(p12[0]));
    EXPECT_EQ(1ul, p12.refs(p12[1]));
    EXPECT_EQ(1ul, p12.refs(p12[2]));
    EXPECT_EQ(1ul, p12.refs(p12[3]));
    EXPECT_EQ(1ul, p12.refs(p12[4]));
    EXPECT_EQ(1ul, p12.refs(p12[6]));
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p12[0]->onedrawsum());
    EXPECT_EQ(3, p12[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p12[1]->onedrawsum());
    EXPECT_EQ(5, p12[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p12[2]->onedrawsum());
    EXPECT_EQ(6, p12[3]->sum());
    EXPECT_EQ(-6, p12[3]->onedrawsum());
    EXPECT_EQ(3, p12[4]->sum());
    EXPECT_EQ(-3, p12[4]->onedrawsum());
    EXPECT_EQ(0, p12[6]->sum());
    EXPECT_EQ(0, p12[6]->onedrawsum());
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    // Test move constructor
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p13 {std::move(e)};
    EXPECT_EQ(2ul, p13.dims());
    EXPECT_EQ(5ul, p13.par_size());
    EXPECT_EQ(6ul, p13.size());
    EXPECT_EQ(1ul, p13.count(0));
    EXPECT_EQ(1ul, p13.count(1));
    EXPECT_EQ(1ul, p13.count(2));
    EXPECT_EQ(1ul, p13.count(3));
    EXPECT_EQ(1ul, p13.count(4));
    EXPECT_EQ(0ul, p13.count(5));
    EXPECT_EQ(0ul, p13.count(6));
    i = 0;
    for (auto it {p13.cbegin()}; it != p13.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p13.refs(p13[0]));
    EXPECT_EQ(1ul, p13.refs(p13[1]));
    EXPECT_EQ(1ul, p13.refs(p13[2]));
    EXPECT_EQ(1ul, p13.refs(p13[3]));
    EXPECT_EQ(1ul, p13.refs(p13[4]));
    EXPECT_EQ(1ul, p13.refs(p13[6]));
    EXPECT_EQ(2, p13[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p13[0]->onedrawsum());
    EXPECT_EQ(3, p13[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p13[1]->onedrawsum());
    EXPECT_EQ(5, p13[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p13[2]->onedrawsum());
    EXPECT_EQ(6, p13[3]->sum());
    EXPECT_EQ(-6, p13[3]->onedrawsum());
    EXPECT_EQ(3, p13[4]->sum());
    EXPECT_EQ(-3, p13[4]->onedrawsum());
    EXPECT_EQ(0, p13[6]->sum());
    EXPECT_EQ(0, p13[6]->onedrawsum());
    // Test move assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> p14(6, v.begin(), v.end());
    p14 = std::move(p12);
    EXPECT_EQ(2ul, p14.dims());
    EXPECT_EQ(5ul, p14.par_size());
    EXPECT_EQ(6ul, p14.size());
    EXPECT_EQ(1ul, p14.count(0));
    EXPECT_EQ(1ul, p14.count(1));
    EXPECT_EQ(1ul, p14.count(2));
    EXPECT_EQ(1ul, p14.count(3));
    EXPECT_EQ(1ul, p14.count(4));
    EXPECT_EQ(0ul, p14.count(5));
    EXPECT_EQ(0ul, p14.count(6));
    i = 0;
    for (auto it {p14.cbegin()}; it != p14.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p14.refs(p14[0]));
    EXPECT_EQ(1ul, p14.refs(p14[1]));
    EXPECT_EQ(1ul, p14.refs(p14[2]));
    EXPECT_EQ(1ul, p14.refs(p14[3]));
    EXPECT_EQ(1ul, p14.refs(p14[4]));
    EXPECT_EQ(1ul, p14.refs(p14[6]));
    EXPECT_EQ(2, p14[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p14[0]->onedrawsum());
    EXPECT_EQ(3, p14[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p14[1]->onedrawsum());
    EXPECT_EQ(5, p14[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p14[2]->onedrawsum());
    EXPECT_EQ(6, p14[3]->sum());
    EXPECT_EQ(-6, p14[3]->onedrawsum());
    EXPECT_EQ(3, p14[4]->sum());
    EXPECT_EQ(-3, p14[4]->onedrawsum());
    EXPECT_EQ(0, p14[6]->sum());
    EXPECT_EQ(0, p14[6]->onedrawsum());
}

TEST(ProjectionTest, ProjectionWorksWithOneDrawAndVanilla) {
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    Wheelwright::Projection<int, Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA, int> p(7);
    //set the seed for reproducibility.
    p.seed(1);
    // We can use the onedrawsum() function of the vectors after each call to
    // p.end_update() to test that the Rademacher averages have been set
    // correctly.
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    EXPECT_EQ(nullptr, p[0]);
    std::size_t i {0};
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test populating the Projection (simple case to then test drop_catchall
    // and clear. A more complex test is done later)
    p.begin_update(0, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(2ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(6ul, p.refs(p[7]));
    EXPECT_NE(p[7]->hash(), p[0]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // Test dims(d)
    p.dims(2);
    EXPECT_EQ(2ul, p.dims());
    // Test drop_catchall()
    p.drop_catchall();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(1ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(nullptr, p[7]);
    // Test clear()
    p.clear();
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    EXPECT_EQ(nullptr, p[0]);
    // Test starting to populate the Projection using the
    // begin_update()-update()-end_update() sequence.
    p.begin_update(0, 1);
    p.update(0, 1);
    p.update(1, 1);
    p.update(2, 1);
    p.update(3, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(0, p[0]->elltwosquared());
    EXPECT_EQ(2, p[3]->sum());
    EXPECT_EQ(2, p[3]->onedrawsum());
    EXPECT_EQ(1, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // update() of parameter that is the only parameter mapped to its vector
    p.begin_update(1, 1);
    p.update(3, 1);
    p.end_update();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0, p[7]->elltwosquared());
    // update() with new parameter causing a new vector to be created.
    p.begin_update(2, 2);
    p.update(4, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(4ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5, p[7]->elltwosquared());
    // Test update() with parameter which before the update is mapped to a
    // vector to which other parameters are mapped, and after the update is
    // mapped to a new vector.
    p.begin_update(4, 2);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(6ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[1], p[3]);
    EXPECT_NE(p[1], p[4]);
    EXPECT_NE(p[1], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[1]->hash(), p[3]->hash());
    EXPECT_NE(p[1]->hash(), p[4]->hash());
    EXPECT_NE(p[1]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[7]->elltwosquared());
    // Test centralized update dropping the cav vector when needed.
    p.begin_update(6, 1);
    p.update(5, 5);
    p.update(6, 5);
    p.end_update();
    EXPECT_EQ(7ul, p.dims());
    EXPECT_EQ(7ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(1ul, p.count(5));
    EXPECT_EQ(1ul, p.count(6));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[5]));
    EXPECT_EQ(2ul, p.refs(p[6]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[5]);
    EXPECT_NE(p[0], p[6]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[2], p[3]);
    EXPECT_NE(p[2], p[4]);
    EXPECT_NE(p[2], p[5]);
    EXPECT_NE(p[2], p[6]);
    EXPECT_NE(p[2], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[3], p[5]);
    EXPECT_NE(p[3], p[6]);
    EXPECT_NE(p[3], p[7]);
    EXPECT_NE(p[4], p[5]);
    EXPECT_NE(p[4], p[6]);
    EXPECT_NE(p[4], p[7]);
    EXPECT_EQ(p[5], p[6]);
    EXPECT_NE(p[5], p[7]);
    EXPECT_EQ(nullptr, p[7]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[5]->hash());
    EXPECT_NE(p[0]->hash(), p[6]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[2]->hash(), p[3]->hash());
    EXPECT_NE(p[2]->hash(), p[4]->hash());
    EXPECT_NE(p[2]->hash(), p[5]->hash());
    EXPECT_NE(p[2]->hash(), p[6]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[3]->hash(), p[5]->hash());
    EXPECT_NE(p[3]->hash(), p[6]->hash());
    EXPECT_NE(p[4]->hash(), p[5]->hash());
    EXPECT_NE(p[4]->hash(), p[6]->hash());
    EXPECT_EQ(p[5]->hash(), p[6]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0 + 1, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0 + 1, p[4]->elltwosquared());
    EXPECT_EQ(5, p[5]->sum());
    EXPECT_DOUBLE_EQ(5, p[5]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 16, p[5]->elltwosquared());
    EXPECT_EQ(5, p[6]->sum());
    EXPECT_DOUBLE_EQ(5, p[6]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 16, p[6]->elltwosquared());
    // Test merge()  -- 1
    // Merge two Projections such that:
    //  * one has the catch-all vector and the other doesn't;
    //  * there is a parameter that is mapped to the catch-all vector in one of
    //    them, and to a non-catch-all vector in the other.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p1(2);
    p1.seed(1);
    p1.begin_update(0, 1);
    p1.update(0, 1);
    p1.end_update();
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(1, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p2(2);
    p2.seed(2);
    p2.begin_update(1, 1);
    p2.update(0, 1);
    p2.update(1, 1);
    p2.end_update();
    EXPECT_EQ(1, p2[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[0]->onedrawsum());
    EXPECT_EQ(1, p2[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[1]->onedrawsum());
    p1.merge(p2);
    EXPECT_EQ(2ul, p1.dims());
    EXPECT_EQ(2ul, p1.par_size());
    EXPECT_EQ(2ul, p1.size());
    EXPECT_EQ(1ul, p1.count(0));
    EXPECT_EQ(1ul, p1.count(1));
    EXPECT_EQ(0ul, p1.count(2));
    i = 0;
    for (auto it {p1.cbegin()}; it != p1.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p1.refs(p1[0]));
    EXPECT_EQ(1ul, p1.refs(p1[1]));
    EXPECT_NE(p1[0], p1[1]);
    EXPECT_EQ(nullptr, p1[2]);
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(2, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p1[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.25, p1[0]->elltwosquared());
    EXPECT_EQ(1, p1[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.25, p1[1]->elltwosquared());
    // Test merge() -- 2
    // Symmetric of the previous test.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p3(2);
    p3.seed(2);
    p3.begin_update(1, 1);
    p3.update(0, 1);
    p3.update(1, 1);
    p3.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p4(2);
    p4.seed(1);
    p4.begin_update(0, 1);
    p4.update(0, 1);
    p4.end_update();
    p3.merge(p4);
    EXPECT_EQ(2ul, p3.dims());
    EXPECT_EQ(2ul, p3.par_size());
    EXPECT_EQ(2ul, p3.size());
    EXPECT_EQ(1ul, p3.count(0));
    EXPECT_EQ(1ul, p3.count(1));
    EXPECT_EQ(0ul, p3.count(2));
    i = 0;
    for (auto it {p3.cbegin()}; it != p3.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p3.refs(p3[0]));
    EXPECT_EQ(1ul, p3.refs(p3[1]));
    EXPECT_NE(p3[0], p3[1]);
    EXPECT_EQ(nullptr, p3[2]);
    EXPECT_NE(p3[0]->hash(), p3[1]->hash());
    EXPECT_EQ(2, p3[0]->sum());
    EXPECT_EQ(-2, p3[0]->onedrawsum());
    EXPECT_EQ(0.25, p3[0]->elltwosquared());
    EXPECT_EQ(1, p3[1]->sum());
    EXPECT_EQ(-1, p3[1]->onedrawsum());
    EXPECT_EQ(0.25, p3[1]->elltwosquared());
    // Test merge() -- 3
    // Merge two Projections both containing the catch-all vector but such that
    // the merged one does not contain the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p5(2);
    p5.seed(1);
    p5.begin_update(0, 1);
    p5.update(0, 1);
    p5.end_update();
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p6(2);
    p6.seed(2);
    p6.begin_update(1, 1);
    p6.update(1, 1);
    p6.end_update();
    EXPECT_EQ(1, p6[1]->sum());
    EXPECT_EQ(-1, p6[1]->onedrawsum());
    p5.merge(p6);
    EXPECT_EQ(2ul, p5.dims());
    EXPECT_EQ(2ul, p5.par_size());
    EXPECT_EQ(2ul, p5.size());
    EXPECT_EQ(1ul, p5.count(0));
    EXPECT_EQ(1ul, p5.count(1));
    EXPECT_EQ(0ul, p5.count(2));
    i = 0;
    for (auto it {p5.cbegin()}; it != p5.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p5.refs(p5[0]));
    EXPECT_EQ(1ul, p5.refs(p5[1]));
    EXPECT_NE(p5[0], p5[1]);
    EXPECT_EQ(nullptr, p5[2]);
    EXPECT_NE(p5[0]->hash(), p5[1]->hash());
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    EXPECT_EQ(0.5, p5[0]->elltwosquared());
    EXPECT_EQ(1, p5[1]->sum());
    EXPECT_EQ(-1, p5[1]->onedrawsum());
    EXPECT_EQ(0.5, p5[1]->elltwosquared());
    // Test merge() -- 4
    // Merge two Projections both containing the catch-all vector and such that
    // the merged one also contains the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p7(2);
    p7.seed(1);
    p7.begin_update(0, 1);
    p7.update(0, 1);
    p7.end_update();
    EXPECT_EQ(1, p7[0]->sum());
    EXPECT_EQ(-1, p7[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p8(2);
    p8.seed(2);
    p8.begin_update(1, 1);
    p8.update(0, 1);
    p8.end_update();
    EXPECT_EQ(1, p8[0]->sum());
    EXPECT_EQ(-1, p8[0]->onedrawsum());
    p7.merge(p8);
    EXPECT_EQ(2ul, p7.dims());
    EXPECT_EQ(1ul, p7.par_size());
    EXPECT_EQ(2ul, p7.size());
    EXPECT_EQ(1ul, p7.count(0));
    EXPECT_EQ(0ul, p7.count(1));
    EXPECT_EQ(0ul, p7.count(2));
    i = 0;
    for (auto it {p7.cbegin()}; it != p7.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p7.refs(p7[0]));
    EXPECT_EQ(1ul, p7.refs(p7[2]));
    EXPECT_NE(p7[0], p7[2]);
    EXPECT_EQ(2, p7[0]->sum());
    EXPECT_EQ(-2, p7[0]->onedrawsum());
    EXPECT_EQ(0.5, p7[0]->elltwosquared());
    EXPECT_EQ(0, p7[2]->sum());
    EXPECT_EQ(0, p7[2]->onedrawsum());
    EXPECT_EQ(0.5, p7[2]->elltwosquared());
    // Test merge()  -- 5
    // Merge two Projections such that a pair of parameters is mapped to the
    // same vector in one Projection and to two different vectors in the other.
    // Vice versa for another pair of parameters.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p9(6);
    p9.seed(1);
    p9.begin_update(0, 1);
    p9.update(0, 1);
    p9.update(1, 1);
    p9.update(2, 2);
    p9.update(3, 3);
    p9.update(4, 3);
    p9.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p10(6);
    p10.seed(2);
    p10.begin_update(1, 1);
    p10.update(0, 1);
    p10.update(1, 2);
    p10.update(2, 3);
    p10.update(3, 3);
    p10.end_update();
    p9.merge(p10);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge()  -- 6
    // Test the special case that the oth Projection is empty.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> e(6);
    p9.merge(e);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge() -- 7
    // Test the special case that this Projection is empty.
    e.merge(p9);
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy constructor;
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p11 {e};
    EXPECT_EQ(2ul, p11.dims());
    EXPECT_EQ(5ul, p11.par_size());
    EXPECT_EQ(6ul, p11.size());
    EXPECT_EQ(1ul, p11.count(0));
    EXPECT_EQ(1ul, p11.count(1));
    EXPECT_EQ(1ul, p11.count(2));
    EXPECT_EQ(1ul, p11.count(3));
    EXPECT_EQ(1ul, p11.count(4));
    EXPECT_EQ(0ul, p11.count(5));
    EXPECT_EQ(0ul, p11.count(6));
    i = 0;
    for (auto it {p11.cbegin()}; it != p11.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p11.refs(p11[0]));
    EXPECT_EQ(1ul, p11.refs(p11[1]));
    EXPECT_EQ(1ul, p11.refs(p11[2]));
    EXPECT_EQ(1ul, p11.refs(p11[3]));
    EXPECT_EQ(1ul, p11.refs(p11[4]));
    EXPECT_EQ(1ul, p11.refs(p11[6]));
    EXPECT_EQ(2, p11[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p11[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p11[0]->elltwosquared());
    EXPECT_EQ(3, p11[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p11[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p11[1]->elltwosquared());
    EXPECT_EQ(5, p11[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p11[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p11[2]->elltwosquared());
    EXPECT_EQ(6, p11[3]->sum());
    EXPECT_EQ(-6, p11[3]->onedrawsum());
    EXPECT_EQ(4.5, p11[3]->elltwosquared());
    EXPECT_EQ(3, p11[4]->sum());
    EXPECT_EQ(-3, p11[4]->onedrawsum());
    EXPECT_EQ(4.5, p11[4]->elltwosquared());
    EXPECT_EQ(0, p11[6]->sum());
    EXPECT_EQ(0, p11[6]->onedrawsum());
    EXPECT_EQ(4.5, p11[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.par_cbegin()}; it != e.par_cend(); ++it)
        ++i;
    EXPECT_EQ(e.par_size(), i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA>  p12(6);
    p12 = e;
    EXPECT_EQ(2ul, p12.dims());
    EXPECT_EQ(5ul, p12.par_size());
    EXPECT_EQ(6ul, p12.size());
    EXPECT_EQ(1ul, p12.count(0));
    EXPECT_EQ(1ul, p12.count(1));
    EXPECT_EQ(1ul, p12.count(2));
    EXPECT_EQ(1ul, p12.count(3));
    EXPECT_EQ(1ul, p12.count(4));
    EXPECT_EQ(0ul, p12.count(5));
    EXPECT_EQ(0ul, p12.count(6));
    i = 0;
    for (auto it {p12.cbegin()}; it != p12.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p12.refs(p12[0]));
    EXPECT_EQ(1ul, p12.refs(p12[1]));
    EXPECT_EQ(1ul, p12.refs(p12[2]));
    EXPECT_EQ(1ul, p12.refs(p12[3]));
    EXPECT_EQ(1ul, p12.refs(p12[4]));
    EXPECT_EQ(1ul, p12.refs(p12[6]));
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p12[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p12[0]->elltwosquared());
    EXPECT_EQ(3, p12[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p12[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p12[1]->elltwosquared());
    EXPECT_EQ(5, p12[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p12[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p12[2]->elltwosquared());
    EXPECT_EQ(6, p12[3]->sum());
    EXPECT_EQ(-6, p12[3]->onedrawsum());
    EXPECT_EQ(4.5, p12[3]->elltwosquared());
    EXPECT_EQ(3, p12[4]->sum());
    EXPECT_EQ(-3, p12[4]->onedrawsum());
    EXPECT_EQ(4.5, p12[4]->elltwosquared());
    EXPECT_EQ(0, p12[6]->sum());
    EXPECT_EQ(0, p12[6]->onedrawsum());
    EXPECT_EQ(4.5, p12[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test move constructor
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p13 {std::move(e)};
    EXPECT_EQ(2ul, p13.dims());
    EXPECT_EQ(5ul, p13.par_size());
    EXPECT_EQ(6ul, p13.size());
    EXPECT_EQ(1ul, p13.count(0));
    EXPECT_EQ(1ul, p13.count(1));
    EXPECT_EQ(1ul, p13.count(2));
    EXPECT_EQ(1ul, p13.count(3));
    EXPECT_EQ(1ul, p13.count(4));
    EXPECT_EQ(0ul, p13.count(5));
    EXPECT_EQ(0ul, p13.count(6));
    i = 0;
    for (auto it {p13.cbegin()}; it != p13.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p13.refs(p13[0]));
    EXPECT_EQ(1ul, p13.refs(p13[1]));
    EXPECT_EQ(1ul, p13.refs(p13[2]));
    EXPECT_EQ(1ul, p13.refs(p13[3]));
    EXPECT_EQ(1ul, p13.refs(p13[4]));
    EXPECT_EQ(1ul, p13.refs(p13[6]));
    EXPECT_EQ(2, p13[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p13[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p13[0]->elltwosquared());
    EXPECT_EQ(3, p13[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p13[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p13[1]->elltwosquared());
    EXPECT_EQ(5, p13[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p13[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p13[2]->elltwosquared());
    EXPECT_EQ(6, p13[3]->sum());
    EXPECT_EQ(-6, p13[3]->onedrawsum());
    EXPECT_EQ(4.5, p13[3]->elltwosquared());
    EXPECT_EQ(3, p13[4]->sum());
    EXPECT_EQ(-3, p13[4]->onedrawsum());
    EXPECT_EQ(4.5, p13[4]->elltwosquared());
    EXPECT_EQ(0, p13[6]->sum());
    EXPECT_EQ(0, p13[6]->onedrawsum());
    EXPECT_EQ(4.5, p13[6]->elltwosquared());
    // Test move assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA> p14(6);
    p14 = std::move(p12);
    EXPECT_EQ(2ul, p14.dims());
    EXPECT_EQ(5ul, p14.par_size());
    EXPECT_EQ(6ul, p14.size());
    EXPECT_EQ(1ul, p14.count(0));
    EXPECT_EQ(1ul, p14.count(1));
    EXPECT_EQ(1ul, p14.count(2));
    EXPECT_EQ(1ul, p14.count(3));
    EXPECT_EQ(1ul, p14.count(4));
    EXPECT_EQ(0ul, p14.count(5));
    EXPECT_EQ(0ul, p14.count(6));
    i = 0;
    for (auto it {p14.cbegin()}; it != p14.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p14.refs(p14[0]));
    EXPECT_EQ(1ul, p14.refs(p14[1]));
    EXPECT_EQ(1ul, p14.refs(p14[2]));
    EXPECT_EQ(1ul, p14.refs(p14[3]));
    EXPECT_EQ(1ul, p14.refs(p14[4]));
    EXPECT_EQ(1ul, p14.refs(p14[6]));
    EXPECT_EQ(2, p14[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p14[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p14[0]->elltwosquared());
    EXPECT_EQ(3, p14[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p14[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p14[1]->elltwosquared());
    EXPECT_EQ(5, p14[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p14[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p14[2]->elltwosquared());
    EXPECT_EQ(6, p14[3]->sum());
    EXPECT_EQ(-6, p14[3]->onedrawsum());
    EXPECT_EQ(4.5, p14[3]->elltwosquared());
    EXPECT_EQ(3, p14[4]->sum());
    EXPECT_EQ(-3, p14[4]->onedrawsum());
    EXPECT_EQ(4.5, p14[4]->elltwosquared());
    EXPECT_EQ(0, p14[6]->sum());
    EXPECT_EQ(0, p14[6]->onedrawsum());
    EXPECT_EQ(4.5, p14[6]->elltwosquared());
    EXPECT_EQ(2, p14[0]->sum());
}

TEST(ProjectionTest, ProjectionWorksWithOneDrawAndVanillaAndDiscrete) {
    // Many of the test cases also check dims(), par_size(), refs(), size(),
    // count(), and operator[]
    // Constructor
    std::set<int> v {1, 2, 3, 4, 5};
    Wheelwright::Projection<int, Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE, int>
        p(7, v.begin(), v.end());
    //set the seed for reproducibility.
    p.seed(1);
    // We can use the onedrawsum() function of the vectors after each call to
    // p.end_update() to test that the Rademacher averages have been set
    // correctly.
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    EXPECT_EQ(nullptr, p[0]);
    std::size_t i {0};
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    // Test populating the Projection (simple case to then test drop_catchall
    // and clear. A more complex test is done later)
    p.begin_update(0, 1);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(2ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(6ul, p.refs(p[7]));
    EXPECT_NE(p[7]->hash(), p[0]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it) {
        bool found_one {false};
        bool found_two {false};
        bool found_three {false};
        bool found_four {false};
        bool found_five {false};
        for (auto it2 {it->domain_begin()}; it2 != it->domain_end(); it2++) {
            switch (*it2) {
                case 1:
                    EXPECT_FALSE(found_one);
                    found_one = true;
                    break;
                case 2:
                    EXPECT_FALSE(found_two);
                    found_two = true;
                    break;
                case 3:
                    EXPECT_FALSE(found_three);
                    found_three = true;
                    break;
                case 4:
                    EXPECT_FALSE(found_four);
                    found_four = true;
                    break;
                case 5:
                    EXPECT_FALSE(found_five);
                    found_five = true;
                    break;
                default:
                    EXPECT_FALSE(true);
                    break;
            }
        }
        EXPECT_TRUE(found_one);
        EXPECT_TRUE(found_two);
        EXPECT_TRUE(found_three);
        EXPECT_TRUE(found_four);
        EXPECT_TRUE(found_five);
        ++i;
    }
    EXPECT_EQ(2ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // Test dims(d)
    p.dims(2);
    EXPECT_EQ(2ul, p.dims());
    // Test drop_catchall()
    p.drop_catchall();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(1ul, p.par_size());
    EXPECT_EQ(1ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(1ul, p.refs(p[0]));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(1ul, i);
    EXPECT_EQ(2, p[0]->sum());
    EXPECT_EQ(-2, p[0]->onedrawsum());
    EXPECT_EQ(1, p[0]->elltwosquared());
    EXPECT_EQ(nullptr, p[7]);
    // Test clear()
    p.clear();
    EXPECT_EQ(0ul, p.dims());
    EXPECT_EQ(0ul, p.par_size());
    EXPECT_EQ(0ul, p.size());
    EXPECT_EQ(0ul, p.count(0));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(0ul, i);
    EXPECT_EQ(nullptr, p[0]);
    // Test starting to populate the Projection using the
    // begin_update()-update()-end_update() sequence.
    p.begin_update(0, 1);
    p.update(0, 1);
    p.update(1, 1);
    p.update(2, 1);
    p.update(3, 2);
    p.end_update();
    EXPECT_EQ(1ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_EQ(1, p[0]->onedrawsum());
    EXPECT_EQ(0, p[0]->elltwosquared());
    EXPECT_EQ(2, p[3]->sum());
    EXPECT_EQ(2, p[3]->onedrawsum());
    EXPECT_EQ(1, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_EQ(0, p[7]->onedrawsum());
    EXPECT_EQ(1, p[7]->elltwosquared());
    // update() of parameter that is the only parameter mapped to its vector
    p.begin_update(1, 1);
    p.update(3, 1);
    p.end_update();
    EXPECT_EQ(2ul, p.dims());
    EXPECT_EQ(4ul, p.par_size());
    EXPECT_EQ(3ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(3ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(3ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0, p[3]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0, p[7]->elltwosquared());
    // update() with new parameter causing a new vector to be created.
    p.begin_update(2, 2);
    p.update(4, 2);
    p.end_update();
    EXPECT_EQ(4ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(4ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(4ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(3ul, p.refs(p[0]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_EQ(p[0], p[1]);
    EXPECT_EQ(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_EQ(p[0]->hash(), p[1]->hash());
    EXPECT_EQ(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(1, p[0]->sum());
    EXPECT_DOUBLE_EQ(1, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5, p[0]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5, p[7]->elltwosquared());
    // Test update() with parameter which before the update is mapped to a
    // vector to which other parameters are mapped, and after the update is
    // mapped to a new vector.
    p.begin_update(4, 2);
    p.update(0, 2);
    p.end_update();
    EXPECT_EQ(6ul, p.dims());
    EXPECT_EQ(5ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[7]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[1], p[3]);
    EXPECT_NE(p[1], p[4]);
    EXPECT_NE(p[1], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[7], p[3]);
    EXPECT_NE(p[7], p[4]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[7]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[1]->hash(), p[3]->hash());
    EXPECT_NE(p[1]->hash(), p[4]->hash());
    EXPECT_NE(p[1]->hash(), p[7]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[7]->hash(), p[3]->hash());
    EXPECT_NE(p[7]->hash(), p[4]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0, p[4]->elltwosquared());
    EXPECT_EQ(0, p[7]->sum());
    EXPECT_DOUBLE_EQ(0, p[7]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0, p[7]->elltwosquared());
    // Test centralized update dropping the cav vector when needed.
    p.begin_update(6, 1);
    p.update(5, 5);
    p.update(6, 5);
    p.end_update();
    EXPECT_EQ(7ul, p.dims());
    EXPECT_EQ(7ul, p.par_size());
    EXPECT_EQ(5ul, p.size());
    EXPECT_EQ(1ul, p.count(0));
    EXPECT_EQ(1ul, p.count(1));
    EXPECT_EQ(1ul, p.count(2));
    EXPECT_EQ(1ul, p.count(3));
    EXPECT_EQ(1ul, p.count(4));
    EXPECT_EQ(1ul, p.count(5));
    EXPECT_EQ(1ul, p.count(6));
    EXPECT_EQ(0ul, p.count(7));
    i = 0;
    for (auto it {p.cbegin()}; it != p.cend(); ++it)
        ++i;
    EXPECT_EQ(5ul, i);
    i = 0;
    for (auto it {p.par_cbegin()}; it != p.par_cend(); ++it)
        ++i;
    EXPECT_EQ(p.par_size(), i);
    EXPECT_EQ(1ul, p.refs(p[0]));
    EXPECT_EQ(2ul, p.refs(p[2]));
    EXPECT_EQ(1ul, p.refs(p[3]));
    EXPECT_EQ(1ul, p.refs(p[4]));
    EXPECT_EQ(2ul, p.refs(p[5]));
    EXPECT_EQ(2ul, p.refs(p[6]));
    EXPECT_NE(p[0], p[1]);
    EXPECT_NE(p[0], p[2]);
    EXPECT_NE(p[0], p[3]);
    EXPECT_NE(p[0], p[4]);
    EXPECT_NE(p[0], p[5]);
    EXPECT_NE(p[0], p[6]);
    EXPECT_NE(p[0], p[7]);
    EXPECT_EQ(p[1], p[2]);
    EXPECT_NE(p[2], p[3]);
    EXPECT_NE(p[2], p[4]);
    EXPECT_NE(p[2], p[5]);
    EXPECT_NE(p[2], p[6]);
    EXPECT_NE(p[2], p[7]);
    EXPECT_NE(p[3], p[4]);
    EXPECT_NE(p[3], p[5]);
    EXPECT_NE(p[3], p[6]);
    EXPECT_NE(p[3], p[7]);
    EXPECT_NE(p[4], p[5]);
    EXPECT_NE(p[4], p[6]);
    EXPECT_NE(p[4], p[7]);
    EXPECT_EQ(p[5], p[6]);
    EXPECT_NE(p[5], p[7]);
    EXPECT_EQ(nullptr, p[7]);
    EXPECT_NE(p[0]->hash(), p[1]->hash());
    EXPECT_NE(p[0]->hash(), p[2]->hash());
    EXPECT_NE(p[0]->hash(), p[3]->hash());
    EXPECT_NE(p[0]->hash(), p[4]->hash());
    EXPECT_NE(p[0]->hash(), p[5]->hash());
    EXPECT_NE(p[0]->hash(), p[6]->hash());
    EXPECT_EQ(p[1]->hash(), p[2]->hash());
    EXPECT_NE(p[2]->hash(), p[3]->hash());
    EXPECT_NE(p[2]->hash(), p[4]->hash());
    EXPECT_NE(p[2]->hash(), p[5]->hash());
    EXPECT_NE(p[2]->hash(), p[6]->hash());
    EXPECT_NE(p[3]->hash(), p[4]->hash());
    EXPECT_NE(p[3]->hash(), p[5]->hash());
    EXPECT_NE(p[3]->hash(), p[6]->hash());
    EXPECT_NE(p[4]->hash(), p[5]->hash());
    EXPECT_NE(p[4]->hash(), p[6]->hash());
    EXPECT_EQ(p[5]->hash(), p[6]->hash());
    EXPECT_EQ(5, p[0]->sum());
    EXPECT_DOUBLE_EQ(-3, p[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 / 9.0 + 0.5 + 2 * 64.0 / 25.0 + 1, p[0]->elltwosquared());
    EXPECT_EQ(1, p[1]->sum());
    EXPECT_EQ(1, p[1]->onedrawsum());
    EXPECT_EQ(1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[1]->elltwosquared());
    EXPECT_EQ(3, p[3]->sum());
    EXPECT_DOUBLE_EQ(3, p[3]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 4.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 1, p[3]->elltwosquared());
    EXPECT_EQ(4, p[4]->sum());
    EXPECT_DOUBLE_EQ(0, p[4]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 9.0 / 2.0 + 2 * 4.0 / 25.0 + 1, p[4]->elltwosquared());
    EXPECT_EQ(5, p[5]->sum());
    EXPECT_DOUBLE_EQ(5, p[5]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 16, p[5]->elltwosquared());
    EXPECT_EQ(5, p[6]->sum());
    EXPECT_DOUBLE_EQ(5, p[6]->onedrawsum());
    EXPECT_DOUBLE_EQ(1.0 + 1.0 / 9.0 + 0.5 + 2 * 4.0 / 25.0 + 16, p[6]->elltwosquared());
    // Test merge()  -- 1
    // Merge two Projections such that:
    //  * one has the catch-all vector and the other doesn't;
    //  * there is a parameter that is mapped to the catch-all vector in one of
    //    them, and to a non-catch-all vector in the other.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p1(2,
                v.begin(), v.end());
    p1.seed(1);
    p1.begin_update(0, 1);
    p1.update(0, 1);
    p1.end_update();
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(1, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p2(2,
                v.begin(), v.end());
    p2.seed(2);
    p2.begin_update(1, 1);
    p2.update(0, 1);
    p2.update(1, 1);
    p2.end_update();
    EXPECT_EQ(1, p2[0]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[0]->onedrawsum());
    EXPECT_EQ(1, p2[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p2[1]->onedrawsum());
    p1.merge(p2);
    EXPECT_EQ(2ul, p1.dims());
    EXPECT_EQ(2ul, p1.par_size());
    EXPECT_EQ(2ul, p1.size());
    EXPECT_EQ(1ul, p1.count(0));
    EXPECT_EQ(1ul, p1.count(1));
    EXPECT_EQ(0ul, p1.count(2));
    i = 0;
    for (auto it {p1.cbegin()}; it != p1.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p1.refs(p1[0]));
    EXPECT_EQ(1ul, p1.refs(p1[1]));
    EXPECT_NE(p1[0], p1[1]);
    EXPECT_EQ(nullptr, p1[2]);
    EXPECT_NE(p1[0]->hash(), p1[1]->hash());
    EXPECT_EQ(2, p1[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p1[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.25, p1[0]->elltwosquared());
    EXPECT_EQ(1, p1[1]->sum());
    EXPECT_DOUBLE_EQ(-1, p1[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.25, p1[1]->elltwosquared());
    // Test merge() -- 2
    // Symmetric of the previous test.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p3(2,
                v.begin(), v.end());
    p3.seed(2);
    p3.begin_update(1, 1);
    p3.update(0, 1);
    p3.update(1, 1);
    p3.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p4(2,
                v.begin(), v.end());
    p4.seed(1);
    p4.begin_update(0, 1);
    p4.update(0, 1);
    p4.end_update();
    p3.merge(p4);
    EXPECT_EQ(2ul, p3.dims());
    EXPECT_EQ(2ul, p3.par_size());
    EXPECT_EQ(2ul, p3.size());
    EXPECT_EQ(1ul, p3.count(0));
    EXPECT_EQ(1ul, p3.count(1));
    EXPECT_EQ(0ul, p3.count(2));
    i = 0;
    for (auto it {p3.cbegin()}; it != p3.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p3.refs(p3[0]));
    EXPECT_EQ(1ul, p3.refs(p3[1]));
    EXPECT_NE(p3[0], p3[1]);
    EXPECT_EQ(nullptr, p3[2]);
    EXPECT_NE(p3[0]->hash(), p3[1]->hash());
    EXPECT_EQ(2, p3[0]->sum());
    EXPECT_EQ(-2, p3[0]->onedrawsum());
    EXPECT_EQ(0.25, p3[0]->elltwosquared());
    EXPECT_EQ(1, p3[1]->sum());
    EXPECT_EQ(-1, p3[1]->onedrawsum());
    EXPECT_EQ(0.25, p3[1]->elltwosquared());
    // Test merge() -- 3
    // Merge two Projections both containing the catch-all vector but such that
    // the merged one does not contain the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p5(2,
                v.begin(), v.end());
    p5.seed(1);
    p5.begin_update(0, 1);
    p5.update(0, 1);
    p5.end_update();
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p6(2,
                v.begin(), v.end());
    p6.seed(2);
    p6.begin_update(1, 1);
    p6.update(1, 1);
    p6.end_update();
    EXPECT_EQ(1, p6[1]->sum());
    EXPECT_EQ(-1, p6[1]->onedrawsum());
    p5.merge(p6);
    EXPECT_EQ(2ul, p5.dims());
    EXPECT_EQ(2ul, p5.par_size());
    EXPECT_EQ(2ul, p5.size());
    EXPECT_EQ(1ul, p5.count(0));
    EXPECT_EQ(1ul, p5.count(1));
    EXPECT_EQ(0ul, p5.count(2));
    i = 0;
    for (auto it {p5.cbegin()}; it != p5.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p5.refs(p5[0]));
    EXPECT_EQ(1ul, p5.refs(p5[1]));
    EXPECT_NE(p5[0], p5[1]);
    EXPECT_EQ(nullptr, p5[2]);
    EXPECT_NE(p5[0]->hash(), p5[1]->hash());
    EXPECT_EQ(1, p5[0]->sum());
    EXPECT_EQ(-1, p5[0]->onedrawsum());
    EXPECT_EQ(0.5, p5[0]->elltwosquared());
    EXPECT_EQ(1, p5[1]->sum());
    EXPECT_EQ(-1, p5[1]->onedrawsum());
    EXPECT_EQ(0.5, p5[1]->elltwosquared());
    // Test merge() -- 4
    // Merge two Projections both containing the catch-all vector and such that
    // the merged one also contains the catch-all vector.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p7(2,
                v.begin(), v.end());
    p7.seed(1);
    p7.begin_update(0, 1);
    p7.update(0, 1);
    p7.end_update();
    EXPECT_EQ(1, p7[0]->sum());
    EXPECT_EQ(-1, p7[0]->onedrawsum());
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p8(2,
                v.begin(), v.end());
    p8.seed(2);
    p8.begin_update(1, 1);
    p8.update(0, 1);
    p8.end_update();
    EXPECT_EQ(1, p8[0]->sum());
    EXPECT_EQ(-1, p8[0]->onedrawsum());
    p7.merge(p8);
    EXPECT_EQ(2ul, p7.dims());
    EXPECT_EQ(1ul, p7.par_size());
    EXPECT_EQ(2ul, p7.size());
    EXPECT_EQ(1ul, p7.count(0));
    EXPECT_EQ(0ul, p7.count(1));
    EXPECT_EQ(0ul, p7.count(2));
    i = 0;
    for (auto it {p7.cbegin()}; it != p7.cend(); ++it)
        ++i;
    EXPECT_EQ(2ul, i);
    EXPECT_EQ(1ul, p7.refs(p7[0]));
    EXPECT_EQ(1ul, p7.refs(p7[2]));
    EXPECT_NE(p7[0], p7[2]);
    EXPECT_EQ(2, p7[0]->sum());
    EXPECT_EQ(-2, p7[0]->onedrawsum());
    EXPECT_EQ(0.5, p7[0]->elltwosquared());
    EXPECT_EQ(0, p7[2]->sum());
    EXPECT_EQ(0, p7[2]->onedrawsum());
    EXPECT_EQ(0.5, p7[2]->elltwosquared());
    // Test merge()  -- 5
    // Merge two Projections such that a pair of parameters is mapped to the
    // same vector in one Projection and to two different vectors in the other.
    // Vice versa for another pair of parameters.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p9(6,
                v.begin(), v.end());
    p9.seed(1);
    p9.begin_update(0, 1);
    p9.update(0, 1);
    p9.update(1, 1);
    p9.update(2, 2);
    p9.update(3, 3);
    p9.update(4, 3);
    p9.end_update();
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p10(6,
                v.begin(), v.end());
    p10.seed(2);
    p10.begin_update(1, 1);
    p10.update(0, 1);
    p10.update(1, 2);
    p10.update(2, 3);
    p10.update(3, 3);
    p10.end_update();
    p9.merge(p10);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge()  -- 6
    // Test the special case that the oth Projection is empty.
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> e(6,
                v.begin(), v.end());
    p9.merge(e);
    EXPECT_EQ(2ul, p9.dims());
    EXPECT_EQ(5ul, p9.par_size());
    EXPECT_EQ(6ul, p9.size());
    EXPECT_EQ(1ul, p9.count(0));
    EXPECT_EQ(1ul, p9.count(1));
    EXPECT_EQ(1ul, p9.count(2));
    EXPECT_EQ(1ul, p9.count(3));
    EXPECT_EQ(1ul, p9.count(4));
    EXPECT_EQ(0ul, p9.count(5));
    EXPECT_EQ(0ul, p9.count(6));
    i = 0;
    for (auto it {p9.cbegin()}; it != p9.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p9.refs(p9[0]));
    EXPECT_EQ(1ul, p9.refs(p9[1]));
    EXPECT_EQ(1ul, p9.refs(p9[2]));
    EXPECT_EQ(1ul, p9.refs(p9[3]));
    EXPECT_EQ(1ul, p9.refs(p9[4]));
    EXPECT_EQ(1ul, p9.refs(p9[6]));
    EXPECT_EQ(2, p9[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p9[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[0]->elltwosquared());
    EXPECT_EQ(3, p9[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p9[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p9[1]->elltwosquared());
    EXPECT_EQ(5, p9[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p9[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p9[2]->elltwosquared());
    EXPECT_EQ(6, p9[3]->sum());
    EXPECT_EQ(-6, p9[3]->onedrawsum());
    EXPECT_EQ(4.5, p9[3]->elltwosquared());
    EXPECT_EQ(3, p9[4]->sum());
    EXPECT_EQ(-3, p9[4]->onedrawsum());
    EXPECT_EQ(4.5, p9[4]->elltwosquared());
    EXPECT_EQ(0, p9[6]->sum());
    EXPECT_EQ(0, p9[6]->onedrawsum());
    EXPECT_EQ(4.5, p9[6]->elltwosquared());
    // Test merge() -- 7
    // Test the special case that this Projection is empty.
    e.merge(p9);
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy constructor;
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p11 {e};
    EXPECT_EQ(2ul, p11.dims());
    EXPECT_EQ(5ul, p11.par_size());
    EXPECT_EQ(6ul, p11.size());
    EXPECT_EQ(1ul, p11.count(0));
    EXPECT_EQ(1ul, p11.count(1));
    EXPECT_EQ(1ul, p11.count(2));
    EXPECT_EQ(1ul, p11.count(3));
    EXPECT_EQ(1ul, p11.count(4));
    EXPECT_EQ(0ul, p11.count(5));
    EXPECT_EQ(0ul, p11.count(6));
    i = 0;
    for (auto it {p11.cbegin()}; it != p11.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p11.refs(p11[0]));
    EXPECT_EQ(1ul, p11.refs(p11[1]));
    EXPECT_EQ(1ul, p11.refs(p11[2]));
    EXPECT_EQ(1ul, p11.refs(p11[3]));
    EXPECT_EQ(1ul, p11.refs(p11[4]));
    EXPECT_EQ(1ul, p11.refs(p11[6]));
    EXPECT_EQ(2, p11[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p11[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p11[0]->elltwosquared());
    EXPECT_EQ(3, p11[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p11[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p11[1]->elltwosquared());
    EXPECT_EQ(5, p11[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p11[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p11[2]->elltwosquared());
    EXPECT_EQ(6, p11[3]->sum());
    EXPECT_EQ(-6, p11[3]->onedrawsum());
    EXPECT_EQ(4.5, p11[3]->elltwosquared());
    EXPECT_EQ(3, p11[4]->sum());
    EXPECT_EQ(-3, p11[4]->onedrawsum());
    EXPECT_EQ(4.5, p11[4]->elltwosquared());
    EXPECT_EQ(0, p11[6]->sum());
    EXPECT_EQ(0, p11[6]->onedrawsum());
    EXPECT_EQ(4.5, p11[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.par_cbegin()}; it != e.par_cend(); ++it)
        ++i;
    EXPECT_EQ(e.par_size(), i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test copy assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE>  p12(6,
                v.begin(), v.end());
    p12 = e;
    EXPECT_EQ(2ul, p12.dims());
    EXPECT_EQ(5ul, p12.par_size());
    EXPECT_EQ(6ul, p12.size());
    EXPECT_EQ(1ul, p12.count(0));
    EXPECT_EQ(1ul, p12.count(1));
    EXPECT_EQ(1ul, p12.count(2));
    EXPECT_EQ(1ul, p12.count(3));
    EXPECT_EQ(1ul, p12.count(4));
    EXPECT_EQ(0ul, p12.count(5));
    EXPECT_EQ(0ul, p12.count(6));
    i = 0;
    for (auto it {p12.cbegin()}; it != p12.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p12.refs(p12[0]));
    EXPECT_EQ(1ul, p12.refs(p12[1]));
    EXPECT_EQ(1ul, p12.refs(p12[2]));
    EXPECT_EQ(1ul, p12.refs(p12[3]));
    EXPECT_EQ(1ul, p12.refs(p12[4]));
    EXPECT_EQ(1ul, p12.refs(p12[6]));
    EXPECT_EQ(2, p12[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p12[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p12[0]->elltwosquared());
    EXPECT_EQ(3, p12[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p12[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p12[1]->elltwosquared());
    EXPECT_EQ(5, p12[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p12[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p12[2]->elltwosquared());
    EXPECT_EQ(6, p12[3]->sum());
    EXPECT_EQ(-6, p12[3]->onedrawsum());
    EXPECT_EQ(4.5, p12[3]->elltwosquared());
    EXPECT_EQ(3, p12[4]->sum());
    EXPECT_EQ(-3, p12[4]->onedrawsum());
    EXPECT_EQ(4.5, p12[4]->elltwosquared());
    EXPECT_EQ(0, p12[6]->sum());
    EXPECT_EQ(0, p12[6]->onedrawsum());
    EXPECT_EQ(4.5, p12[6]->elltwosquared());
    EXPECT_EQ(2ul, e.dims());
    EXPECT_EQ(5ul, e.par_size());
    EXPECT_EQ(6ul, e.size());
    EXPECT_EQ(1ul, e.count(0));
    EXPECT_EQ(1ul, e.count(1));
    EXPECT_EQ(1ul, e.count(2));
    EXPECT_EQ(1ul, e.count(3));
    EXPECT_EQ(1ul, e.count(4));
    EXPECT_EQ(0ul, e.count(5));
    EXPECT_EQ(0ul, e.count(6));
    i = 0;
    for (auto it {e.cbegin()}; it != e.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, e.refs(e[0]));
    EXPECT_EQ(1ul, e.refs(e[1]));
    EXPECT_EQ(1ul, e.refs(e[2]));
    EXPECT_EQ(1ul, e.refs(e[3]));
    EXPECT_EQ(1ul, e.refs(e[4]));
    EXPECT_EQ(1ul, e.refs(e[6]));
    EXPECT_EQ(2, e[0]->sum());
    EXPECT_DOUBLE_EQ(-2, e[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[0]->elltwosquared());
    EXPECT_EQ(3, e[1]->sum());
    EXPECT_DOUBLE_EQ(-3, e[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, e[1]->elltwosquared());
    EXPECT_EQ(5, e[2]->sum());
    EXPECT_DOUBLE_EQ(-5, e[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, e[2]->elltwosquared());
    EXPECT_EQ(6, e[3]->sum());
    EXPECT_EQ(-6, e[3]->onedrawsum());
    EXPECT_EQ(4.5, e[3]->elltwosquared());
    EXPECT_EQ(3, e[4]->sum());
    EXPECT_EQ(-3, e[4]->onedrawsum());
    EXPECT_EQ(4.5, e[4]->elltwosquared());
    EXPECT_EQ(0, e[6]->sum());
    EXPECT_EQ(0, e[6]->onedrawsum());
    EXPECT_EQ(4.5, e[6]->elltwosquared());
    // Test move constructor
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p13 {std::move(e)};
    EXPECT_EQ(2ul, p13.dims());
    EXPECT_EQ(5ul, p13.par_size());
    EXPECT_EQ(6ul, p13.size());
    EXPECT_EQ(1ul, p13.count(0));
    EXPECT_EQ(1ul, p13.count(1));
    EXPECT_EQ(1ul, p13.count(2));
    EXPECT_EQ(1ul, p13.count(3));
    EXPECT_EQ(1ul, p13.count(4));
    EXPECT_EQ(0ul, p13.count(5));
    EXPECT_EQ(0ul, p13.count(6));
    i = 0;
    for (auto it {p13.cbegin()}; it != p13.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p13.refs(p13[0]));
    EXPECT_EQ(1ul, p13.refs(p13[1]));
    EXPECT_EQ(1ul, p13.refs(p13[2]));
    EXPECT_EQ(1ul, p13.refs(p13[3]));
    EXPECT_EQ(1ul, p13.refs(p13[4]));
    EXPECT_EQ(1ul, p13.refs(p13[6]));
    EXPECT_EQ(2, p13[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p13[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p13[0]->elltwosquared());
    EXPECT_EQ(3, p13[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p13[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p13[1]->elltwosquared());
    EXPECT_EQ(5, p13[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p13[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p13[2]->elltwosquared());
    EXPECT_EQ(6, p13[3]->sum());
    EXPECT_EQ(-6, p13[3]->onedrawsum());
    EXPECT_EQ(4.5, p13[3]->elltwosquared());
    EXPECT_EQ(3, p13[4]->sum());
    EXPECT_EQ(-3, p13[4]->onedrawsum());
    EXPECT_EQ(4.5, p13[4]->elltwosquared());
    EXPECT_EQ(0, p13[6]->sum());
    EXPECT_EQ(0, p13[6]->onedrawsum());
    EXPECT_EQ(4.5, p13[6]->elltwosquared());
    // Test move assignment
    Wheelwright::Projection<int,Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> p14(6,
                v.begin(), v.end());
    p14 = std::move(p12);
    EXPECT_EQ(2ul, p14.dims());
    EXPECT_EQ(5ul, p14.par_size());
    EXPECT_EQ(6ul, p14.size());
    EXPECT_EQ(1ul, p14.count(0));
    EXPECT_EQ(1ul, p14.count(1));
    EXPECT_EQ(1ul, p14.count(2));
    EXPECT_EQ(1ul, p14.count(3));
    EXPECT_EQ(1ul, p14.count(4));
    EXPECT_EQ(0ul, p14.count(5));
    EXPECT_EQ(0ul, p14.count(6));
    i = 0;
    for (auto it {p14.cbegin()}; it != p14.cend(); ++it)
        ++i;
    EXPECT_EQ(6ul, i);
    EXPECT_EQ(1ul, p14.refs(p14[0]));
    EXPECT_EQ(1ul, p14.refs(p14[1]));
    EXPECT_EQ(1ul, p14.refs(p14[2]));
    EXPECT_EQ(1ul, p14.refs(p14[3]));
    EXPECT_EQ(1ul, p14.refs(p14[4]));
    EXPECT_EQ(1ul, p14.refs(p14[6]));
    EXPECT_EQ(2, p14[0]->sum());
    EXPECT_DOUBLE_EQ(-2, p14[0]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p14[0]->elltwosquared());
    EXPECT_EQ(3, p14[1]->sum());
    EXPECT_DOUBLE_EQ(-3, p14[1]->onedrawsum());
    EXPECT_DOUBLE_EQ(0.5, p14[1]->elltwosquared());
    EXPECT_EQ(5, p14[2]->sum());
    EXPECT_DOUBLE_EQ(-5, p14[2]->onedrawsum());
    EXPECT_DOUBLE_EQ(2.5, p14[2]->elltwosquared());
    EXPECT_EQ(6, p14[3]->sum());
    EXPECT_EQ(-6, p14[3]->onedrawsum());
    EXPECT_EQ(4.5, p14[3]->elltwosquared());
    EXPECT_EQ(3, p14[4]->sum());
    EXPECT_EQ(-3, p14[4]->onedrawsum());
    EXPECT_EQ(4.5, p14[4]->elltwosquared());
    EXPECT_EQ(0, p14[6]->sum());
    EXPECT_EQ(0, p14[6]->onedrawsum());
    EXPECT_EQ(4.5, p14[6]->elltwosquared());
    EXPECT_EQ(2, p14[0]->sum());
}
