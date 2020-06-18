/** @file boundstests.cpp
 * Tests for the functions in Bounds.h and Bounds_impl.h
 *
 * @author Matteo Riondato
 * @date 2017 07 31
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
#include <iostream>
#include <list>
#include <vector>

#include "Bounds.h"
#include "Projection.h"

#include "gtest/gtest.h"

TEST(BoundsTest, BoundsWorks) {
    // Test elltwosquared (with arithmetic type)
    EXPECT_DOUBLE_EQ(2, Wheelwright::impl::elltwosquared(2));
    // Test elltwosquared (with Vector type);
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> vl
    {std::sqrt(2), 0, 1};
    EXPECT_DOUBLE_EQ(2, Wheelwright::impl::elltwosquared(vl));
    // Test maxelltwosquared
    std::vector<double> vml {1, 2, 3, 6, 5, 4};
    EXPECT_DOUBLE_EQ(6, Wheelwright::impl::maxelltwosquared(vml.begin(),
                vml.end()));
    std::list<double> lml {vml.begin(), vml.end()};
    EXPECT_DOUBLE_EQ(6, Wheelwright::impl::maxelltwosquared(lml.begin(),
                lml.end()));
    // Test onedrawsum (with arithmetic type)
    EXPECT_DOUBLE_EQ(2, Wheelwright::impl::onedrawsum(2));
    // Test onedrawsum (with Vector type)
    Wheelwright::impl::Vector<double,Wheelwright::Features::ONEDRAW> vo
    {2, 0, 1, std::vector<Wheelwright::impl::Rademacher>(1,
            Wheelwright::impl::Rademacher::NEG)};
    EXPECT_DOUBLE_EQ(-2, Wheelwright::impl::onedrawsum(vo));
    vo.insert(2, 1, 1, std::vector<Wheelwright::impl::Rademacher>(1,
            Wheelwright::impl::Rademacher::POS));
    EXPECT_DOUBLE_EQ(0, Wheelwright::impl::onedrawsum(vo));
    // Test era_bound* functions
    // Test era_bound_MASSART (with values)
    const double val {2 * std::sqrt(std::log(20)) / 10};
    EXPECT_DOUBLE_EQ(val, Wheelwright::era_bound_MASSART(std::sqrt(2), 20, 10));
    // Test era_bound_MASSART (with random-access-iterator)
    std::vector<double> v(19, 1);
    v.push_back(2);
    EXPECT_DOUBLE_EQ(val,
            Wheelwright::era_bound_MASSART(v.begin(), v.end(), 10));
    // Test era_bound_MASSART (without random-access-iterator)
    std::list<double> l {v.begin(), v.end()};
    EXPECT_DOUBLE_EQ(val,
            Wheelwright::era_bound_MASSART(l.begin(), l.end(), 10));
    // Test era_bound_RU (with random-access-iterator for arithmetic type)
    v.clear();
    for (std::size_t i {0}; i < 20; ++i)
        v.push_back(2);
    // The actual bound should be std::sqrt(std::log(20)) / 5
    // = 0.3461636765204571, but it doesn't seem attainable by NLopt, so we use
    // what we can actually get.
    EXPECT_DOUBLE_EQ(0.34616367652049429,
            Wheelwright::era_bound_RU(v.begin(), v.end(), 10, 1e-10));
    // Test era_bound_AGOR (with random-access-iterator for arithmetic type)
    // In this case the bound should be the same as the one for era_bound_RU,
    // but due to the different order of operations, it isn't.
    EXPECT_DOUBLE_EQ(0.34616367652045782,
            Wheelwright::era_bound_AGOR(v.begin(), v.end(), 10, 1e-10));
    // Test era_bound_RU (without random-access-iterator for arithmetic type)
    l.clear();
    l.insert(l.end(), v.begin(), v.end());
    EXPECT_DOUBLE_EQ(0.34616367652049429,
            Wheelwright::era_bound_RU(l.begin(), l.end(), 10, 1e-10));
    // Test era_bound_AGOR (without random-access-iterator for arithmetic type)
    EXPECT_DOUBLE_EQ(0.34616367652045782,
            Wheelwright::era_bound_AGOR(l.begin(), l.end(), 10, 1e-10));
    // Test era_bound* with different values, just to be sure.
    for (std::size_t i {10}; i < 20; ++i)
        v[i] = 1;
    EXPECT_DOUBLE_EQ(0.3177624141418931,
            Wheelwright::era_bound_RU(v.begin(), v.end(), 10, 1e-10));
    EXPECT_DOUBLE_EQ(0.31776241414189299,
            Wheelwright::era_bound_AGOR(v.begin(), v.end(), 10, 1e-10));
    // Test era_bound* functions with Vectors
    std::vector<double> domain {1.0, 2.0};
    std::vector<Wheelwright::impl::Rademacher> radeone {
        Wheelwright::impl::Rademacher::POS};
    std::vector<Wheelwright::impl::Rademacher> radetwo {
        Wheelwright::impl::Rademacher::POS, Wheelwright::impl::Rademacher::POS};
    std::vector<Wheelwright::impl::Rademacher> radethree {
        Wheelwright::impl::Rademacher::POS, Wheelwright::impl::Rademacher::POS,
        Wheelwright::impl::Rademacher::POS};
    std::vector<Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA
        | Wheelwright::Features::DISCRETE | Wheelwright::Features::ONEDRAW>>
        vd;
    vd.emplace_back(domain.cbegin(), domain.cend(), 0, 0, 3, radethree, -2);
    vd.emplace_back(domain.cbegin(), domain.cend(), 1, 0, 1, radeone, -1);
    vd.emplace_back(domain.cbegin(), domain.cend(), 1, 1, 1, radeone, 0);
    vd.emplace_back(domain.cbegin(), domain.cend(), 1, 0, 2, radetwo, 0.5);
    vd.emplace_back(domain.cbegin(), domain.cend(), 1, 0, 3, radethree, 0.7);
    // Test odara() (with random-access-iterator for Vector type)
    EXPECT_DOUBLE_EQ(0.5, Wheelwright::odara(vd.begin(), vd.end(), 3));
    // Test era_bound_AGOR (with random-access-iterator for Vector type)
    EXPECT_DOUBLE_EQ(1.4310440991628803,
            Wheelwright::era_bound_AGOR(vd.begin(), vd.end(), 3, 1e-10));
    // Test era_bound_RU (with random-access-iterator for Vector type)
    EXPECT_DOUBLE_EQ(1.4310440991635938,
            Wheelwright::era_bound_RU(vd.begin(), vd.end(), 3, 1e-10));
    // Test era_bound_AGOR_FD (with random-access-iterator for Vector type)
    EXPECT_DOUBLE_EQ(0.70796233581524881,
            Wheelwright::era_bound_AGOR_FD(vd.begin(), vd.end(), 3,
                domain.begin(), domain.end()));
    // Test era_bound_RU_FD (with random-access-iterator for Vector type)
    EXPECT_DOUBLE_EQ(0.70796233595447067,
            Wheelwright::era_bound_RU_FD(vd.begin(), vd.end(), 3,
                domain.begin(), domain.end()));
    // Test era_bound() (with random-access-iterator for Vector type, but since
    // the pointed-to-type of the iterator or whether it has random access or
    // not does not matter for era_bound(), the following are the only tests for
    // this function)
    EXPECT_DOUBLE_EQ(1.4310440991628803,
            Wheelwright::era_bound(vd.begin(), vd.end(), 3,
                Wheelwright::EraUBMethod::AGOR, 1e-10));
    EXPECT_DOUBLE_EQ(1.4310440991635938,
            Wheelwright::era_bound(vd.begin(), vd.end(), 3,
                Wheelwright::EraUBMethod::RU, 1e-10));
    EXPECT_DOUBLE_EQ(0.70796233581524881,
            Wheelwright::era_bound(vd.begin(), vd.end(), 3,
                Wheelwright::EraUBMethod::AGOR_FD, domain.cbegin(),
                domain.cend(), 1e-10));
    EXPECT_DOUBLE_EQ(0.70796233581524881,
            Wheelwright::era_bound(vd.begin(), vd.end(), 3,
                Wheelwright::EraUBMethod::AGOR_FD, domain.begin(), domain.end(),
                1e-10));
    //EXPECT_DOUBLE_EQ(0.70796233581524948,
    EXPECT_DOUBLE_EQ(0.70796233581524937,
            Wheelwright::era_bound(vd.begin(), vd.end(), 3,
                Wheelwright::EraUBMethod::RU_FD, domain.begin(), domain.end(),
                1e-10));
    EXPECT_DOUBLE_EQ(0.70796233581524937,
            Wheelwright::era_bound(vd.begin(), vd.end(), 3,
                Wheelwright::EraUBMethod::RU_FD, domain.begin(), domain.end(), 1e-10));
    // Test era_bound* with Projection with one combination of template
    // parameters. These tests also act as tests of the era_bound_* functions
    // with non-random-access iterators of Vectors.
    // This Projection does not contain exactly the Vectors that we defined
    // earlier because it is just too much work to create Vectors that have
    // different entries.
    Wheelwright::Projection<int,Wheelwright::Features::VANILLA|
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE> p(6,
                domain.begin(), domain.end());
    p.seed(2);
    p.begin_update(0, 1);
    p.update(1, 1);
    p.update(3, 1);
    p.update(4, 1);
    p.end_update();
    p.begin_update(1, 1);
    p.update(2, 1);
    p.update(3, 1);
    p.update(4, 1);
    p.end_update();
    p.begin_update(2, 1);
    p.update(2, 1);
    p.update(4, 1);
    p.end_update();
    // Test odara()
    EXPECT_DOUBLE_EQ(0.5, Wheelwright::odara(p));
    EXPECT_DOUBLE_EQ(0.51804654246098969,
            Wheelwright::era_bound_AGOR(p.cbegin(), p.cend(), p.dims()));
    EXPECT_DOUBLE_EQ(0.554600277415751,
            Wheelwright::era_bound_MASSART(p.cbegin(), p.cend(), p.dims()));
    EXPECT_DOUBLE_EQ(0.51804654246145809,
            Wheelwright::era_bound_RU(p.cbegin(), p.cend(), p.dims()));
    EXPECT_DOUBLE_EQ(0.73452843186271632,
            Wheelwright::era_bound_AGOR_FD(p.cbegin(), p.cend(), p.dims(),
                p.cbegin()->domain_begin(), p.cbegin()->domain_end()));
    EXPECT_DOUBLE_EQ(0.73452843186270045,
            Wheelwright::era_bound_RU_FD(p.cbegin(), p.cend(), p.dims(),
                p.cbegin()->domain_begin(), p.cbegin()->domain_end()));
    // Test ra_bound* functions. We use the ERA bound from the l.h.s. of the
    // above call to era_bound_RU, and the odara for p.
    // Test impl::ra_bound_MCDIARMIDONEDRAW
    EXPECT_DOUBLE_EQ(1.1375335739364312,
            Wheelwright::impl::ra_bound_MCDIARMIDONEDRAW(0.51804654246145809,
                p.dims(), 0.1, 1));
    EXPECT_DOUBLE_EQ(1.1194870314749732,
            Wheelwright::impl::ra_bound_MCDIARMIDONEDRAW(0.5, p.dims(),
                0.1, 1));
    // Test impl::ra_bound_SELFBOUNDING
    EXPECT_DOUBLE_EQ(1.639977618287531,
            Wheelwright::impl::ra_bound_SELFBOUNDING(0.51804654246145809,
                p.dims(), 0.1, 1));
    // Test ra_bound (with value)
    EXPECT_DOUBLE_EQ(1.1375335739364312,
            Wheelwright::ra_bound(0.51804654246145809, p.dims(), 0.1, 1,
                Wheelwright::RaUBMethod::MCDIARMID));
    EXPECT_DOUBLE_EQ(1.1194870314749732,
            Wheelwright::ra_bound(0.5, p.dims(), 0.1, 1,
                Wheelwright::RaUBMethod::ONEDRAW));
    EXPECT_DOUBLE_EQ(1.639977618287531,
            Wheelwright::ra_bound(0.51804654246145809, p.dims(), 0.1, 1,
                Wheelwright::RaUBMethod::SELFBOUNDING));
    // Test ra_bound (with Projection, thus with non-random-access iterator of
    // Vectors. We do not test with random-access iterator of Vectors or with
    // iterators of numbers because these make a difference only in era_bound.
    // For the same reason we only test with a single value of
    // Wheelwright::EraUBMethod)
    EXPECT_DOUBLE_EQ(1.1375335739364312, Wheelwright::ra_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::MCDIARMID,
                Wheelwright::EraUBMethod::RU));
    EXPECT_DOUBLE_EQ(1.1194870314749732, Wheelwright::ra_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::ONEDRAW,
                Wheelwright::EraUBMethod::RU));
    EXPECT_DOUBLE_EQ(1.639977618287531, Wheelwright::ra_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::SELFBOUNDING,
                Wheelwright::EraUBMethod::RU));
    // Test dev_bound functions
    // Test dev_bound (with value)
    // We use the RA bound from the l.h.s. of the above call to ra_bound with
    // MCDIARMID and RU.
    EXPECT_DOUBLE_EQ(2.98167079367367377,
            Wheelwright::dev_bound(1.1375335739364312, p.dims(), 0.1, 1,
                Wheelwright::MaxDevType::ABSOL));
    EXPECT_DOUBLE_EQ(2.8945541793478354744920,
            Wheelwright::dev_bound(1.1375335739364312, p.dims(), 0.1, 1,
                Wheelwright::MaxDevType::PLAIN));
    // Test dev_bound (with Projection, thus with non-random-access iterator of
    // Vectors. We do not test with random-access iterator of Vectors or with
    // iterators of numbers because these make a difference only in era_bound.
    // We test with different values for RaUBMethod and MaxDevType because they
    // have an impact on dev_bound)
    //
    // It is actually correct that the l.h.s. of the following call is not the
    // same as the l.h.s. of the above call to dev_bound with ABSOL, because the
    // following call uses a single application of the McDiarmid inequality and
    // returns a correct result (i.e., a real bound to the maximum deviation)
    // with probability at least 0.9, while the one above uses two applications
    // of McDiarmid inequality and returns a correct result with probability at
    // least 0.8
    EXPECT_DOUBLE_EQ(3.1559040223253501, Wheelwright::dev_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::MCDIARMID,
                Wheelwright::EraUBMethod::RU,
                Wheelwright::MaxDevType::ABSOL));
    EXPECT_DOUBLE_EQ(2.8945541793478354744920, Wheelwright::dev_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::MCDIARMID,
                Wheelwright::EraUBMethod::RU,
                Wheelwright::MaxDevType::PLAIN));
    EXPECT_DOUBLE_EQ(4.5698833461610091, Wheelwright::dev_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::SELFBOUNDING,
                Wheelwright::EraUBMethod::RU,
                Wheelwright::MaxDevType::ABSOL));
    EXPECT_DOUBLE_EQ(4.4923867162621347, Wheelwright::dev_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::SELFBOUNDING,
                Wheelwright::EraUBMethod::RU,
                Wheelwright::MaxDevType::PLAIN));
    EXPECT_DOUBLE_EQ(2.8584610944249196, Wheelwright::dev_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::ONEDRAW,
                Wheelwright::EraUBMethod::RU,
                Wheelwright::MaxDevType::PLAIN));
    EXPECT_DOUBLE_EQ(3.1198109374024341, Wheelwright::dev_bound(p, 0.1, 1,
                Wheelwright::RaUBMethod::ONEDRAW,
                Wheelwright::EraUBMethod::RU,
                Wheelwright::MaxDevType::ABSOL));
}
