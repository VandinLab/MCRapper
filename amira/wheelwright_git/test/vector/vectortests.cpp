/** @file projvectests.cpp
 * Tests for the Vector class
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

#include <cmath>
#include <cstddef>

#include <boost/functional/hash.hpp>

#include "Projection.h"
#include "Vector.h"

#include "gtest/gtest.h"

TEST(VectorTest, VectorVANILLAWorks) {
    // Test constructor with values (and sum(), elltwosquared(), hash(), and
    // elltwo())
    std::size_t hash {0};
    boost::hash_combine(hash, 0u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 1.0);
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> v(1, 0, 1);
    EXPECT_EQ(1, v.sum());
    EXPECT_EQ(1, v.elltwosquared());
    EXPECT_EQ(1, v.elltwo());
    EXPECT_EQ(hash, v.hash());
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> h(1, 0, 1,
            2);
    EXPECT_EQ(1, h.sum());
    EXPECT_EQ(4, h.elltwosquared());
    EXPECT_EQ(2, h.elltwo());
    EXPECT_EQ(v.hash(), h.hash());
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> h2(1, 0, 1,
            std::numeric_limits<double>::max());
    EXPECT_EQ(1, h2.sum());
    EXPECT_EQ(1, h2.elltwosquared());
    EXPECT_EQ(1, h2.elltwo());
    EXPECT_EQ(v.hash(), h.hash());
    // Test insert()
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2u);
    boost::hash_combine(hash, 2.0);
    v.insert(2, 1, 2);
    EXPECT_EQ(5, v.sum());
    EXPECT_EQ(9, v.elltwosquared());
    EXPECT_EQ(3, v.elltwo());
    EXPECT_EQ(hash, v.hash());
    h.insert(2, 1, 2, 1);
    EXPECT_EQ(5, h.sum());
    EXPECT_EQ(6, h.elltwosquared());
    EXPECT_EQ(std::sqrt(6), h.elltwo());
    h2.insert(2, 1, 2, std::numeric_limits<double>::max());
    EXPECT_EQ(5, h2.sum());
    EXPECT_EQ(9, h2.elltwosquared());
    EXPECT_EQ(3, h2.elltwo());
    // Test constructor with vector and new components
    boost::hash_combine(hash, 3u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2.0);
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> w(v, 2, 3,
            1);
    EXPECT_EQ(7, w.sum());
    EXPECT_EQ(13, w.elltwosquared());
    EXPECT_EQ(std::sqrt(13), w.elltwo());
    EXPECT_EQ(hash, w.hash());
    // Test "merge" constructor
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> z(2, 3, 1);
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> x(v, z);
    EXPECT_EQ(w.sum(), x.sum());
    EXPECT_EQ(w.elltwosquared(), x.elltwosquared());
    EXPECT_EQ(w.elltwo(), x.elltwo());
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> x2(z, v);
    EXPECT_EQ(x.sum(), x2.sum());
    EXPECT_EQ(x.elltwosquared(), x2.elltwosquared());
    EXPECT_EQ(x.elltwo(), x2.elltwo());
    EXPECT_EQ(x.hash(), x2.hash());
    // Test operator<
    EXPECT_FALSE(v < v);
    // Test operator==
    EXPECT_TRUE(v == v);
    EXPECT_FALSE(v == w);
    // Test that the computation of the hash satisfies basic requirements for
    // vectors that are prefixes/suffixes of others.
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> a(1, 0, 1);
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> b(1, 1, 1);
    Wheelwright::impl::Vector<double,Wheelwright::Features::VANILLA> c(1, 0, 2);
    EXPECT_FALSE(a.hash() == b.hash());
    EXPECT_FALSE(a.hash() == c.hash());
    EXPECT_FALSE(b.hash() == c.hash());
}

TEST(VectorTest, VectorDISCRETEWorks) {
    // Test constructor with values (and sum() and hash())
    std::size_t hash {0};
    std::set<int> vals {0, 1, 2, 3};
    boost::hash_combine(hash, 0u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 1);
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> v(
            vals.begin(), vals.end(), 1, 0, 1);
    EXPECT_EQ(1, v.sum());
    EXPECT_EQ(hash, v.hash());
    // Test domain_begin() and domain_end()
    std::size_t i {0};
    bool found_one {false};
    bool found_three {false};
    bool found_two {false};
    for (auto it {v.domain_begin()}; it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test insert()
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2u);
    boost::hash_combine(hash, 2);
    v.insert(2, 1, 2);
    EXPECT_EQ(5, v.sum());
    EXPECT_EQ(hash, v.hash());
    // Test counts_begin() and counts_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.counts_begin(); it != v.counts_end(); ++it) {
        switch(*it) {
            case 0:
                EXPECT_FALSE(found_three);
                found_three = true;
                break;
            case 1:
                EXPECT_FALSE(found_one);
                found_one = true;
                break;
            case 2:
                EXPECT_FALSE(found_two);
                found_two = true;
                break;
            default:
                EXPECT_FALSE(true);
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test domain_begin() and domain_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.domain_begin(); it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test constructor with vector and new components
    boost::hash_combine(hash, 3u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2);
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> w(v, 2, 3, 1);
    EXPECT_EQ(7, w.sum());
    EXPECT_EQ(hash, w.hash());
    // Test "merge" constructor
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> z(vals.begin(), vals.end(), 2, 3, 1);
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> x(v, z);
    EXPECT_EQ(w.sum(), x.sum());
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> x2(z, v);
    EXPECT_EQ(x.sum(), x2.sum());
    EXPECT_EQ(x.hash(), x2.hash());
    // Test operator<
    EXPECT_FALSE(v < v);
    // Test operator==
    EXPECT_TRUE(v == v);
    EXPECT_FALSE(v == w);
    // Test that the computation of the hash satisfies basic requirements for
    // vectors that are prefixes/suffixes of others.
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> a(
            vals.begin(), vals.end(), 1, 0, 1);
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> b(
            vals.begin(), vals.end(), 1, 1, 1);
    Wheelwright::impl::Vector<int,Wheelwright::Features::DISCRETE> c(
            vals.begin(), vals.end(), 1, 0, 2);
    EXPECT_FALSE(a.hash() == b.hash());
    EXPECT_FALSE(a.hash() == c.hash());
    EXPECT_FALSE(b.hash() == c.hash());
}

TEST(VectorTest, VectorONEDRAWWorks) {
    // Test constructor with values (and sum() and hash())
    std::size_t hash {0};
    boost::hash_combine(hash, 0u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 1);
    std::vector<Wheelwright::impl::Rademacher> pos(1,
            Wheelwright::impl::Rademacher::POS);
    std::vector<Wheelwright::impl::Rademacher> neg(1,
            Wheelwright::impl::Rademacher::NEG);
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> v(1, 0, 1,
           neg);
    EXPECT_EQ(1, v.sum());
    EXPECT_EQ(-1, v.onedrawsum());
    EXPECT_EQ(hash, v.hash());
    // Test insert()
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2u);
    boost::hash_combine(hash, 2);
    std::vector<Wheelwright::impl::Rademacher> pos2(2,
            Wheelwright::impl::Rademacher::POS);
    v.insert(2, 1, 2, pos2);
    EXPECT_EQ(5, v.sum());
    EXPECT_EQ(3, v.onedrawsum());
    EXPECT_EQ(hash, v.hash());
    // Test constructor with vector and new components
    boost::hash_combine(hash, 3u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2);
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> w(v, 2, 3, 1,
            neg);
    EXPECT_EQ(7, w.sum());
    EXPECT_EQ(1, w.onedrawsum());
    EXPECT_EQ(hash, w.hash());
    // Test "merge" constructor
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> z(2, 3, 1,
            neg);
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> x(v, z);
    EXPECT_EQ(w.sum(), x.sum());
    EXPECT_EQ(w.onedrawsum(), x.onedrawsum());
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> x2(z, v);
    EXPECT_EQ(x.sum(), x2.sum());
    EXPECT_EQ(x.onedrawsum(), x2.onedrawsum());
    EXPECT_EQ(x.hash(), x2.hash());
    // Test operator<
    EXPECT_FALSE(v < v);
    // Test operator==
    EXPECT_TRUE(v == v);
    EXPECT_FALSE(v == w);
    // Test that the computation of the hash satisfies basic requirements for
    // vectors that are prefixes/suffixes of others.
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> a(1, 0, 1,
            pos);
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> b(1, 1, 1,
            pos);
    Wheelwright::impl::Vector<int,Wheelwright::Features::ONEDRAW> c(1, 0, 2,
            pos2);
    EXPECT_FALSE(a.hash() == b.hash());
    EXPECT_FALSE(a.hash() == c.hash());
    EXPECT_FALSE(b.hash() == c.hash());
}

TEST(VectorTest, VectorVANILLADISCRETEWorks) {
    // Test constructor with values (and sum(), elltwosquared(), hash(), and
    // elltwo())
    std::size_t hash {0};
    std::set<int> vals {0, 1, 2, 3};
    boost::hash_combine(hash, 0u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 1);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> v(
                vals.begin(), vals.end(), 1, 0, 1, 1);
    EXPECT_EQ(1, v.sum());
    EXPECT_EQ(1, v.elltwosquared());
    EXPECT_EQ(1, v.elltwo());
    EXPECT_EQ(hash, v.hash());
    // Test domain_begin() and domain_end()
    std::size_t i {0};
    bool found_one {false};
    bool found_three {false};
    bool found_two {false};
    for (auto it {v.domain_begin()}; it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test insert()
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2u);
    boost::hash_combine(hash, 2);
    v.insert(2, 1, 2, 1);
    EXPECT_EQ(5, v.sum());
    EXPECT_EQ(3, v.elltwosquared());
    EXPECT_EQ(std::sqrt(3), v.elltwo());
    EXPECT_EQ(hash, v.hash());
    // Test counts_begin() and counts_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.counts_begin(); it != v.counts_end(); ++it) {
        switch(*it) {
            case 0:
                EXPECT_FALSE(found_three);
                found_three = true;
                break;
            case 1:
                EXPECT_FALSE(found_one);
                found_one = true;
                break;
            case 2:
                EXPECT_FALSE(found_two);
                found_two = true;
                break;
            default:
                EXPECT_FALSE(true);
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test domain_begin() and domain_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.domain_begin(); it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test constructor with vector and new components
    boost::hash_combine(hash, 3u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2);
    Wheelwright::impl::Vector<int,Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> w(v, 2, 3, 1, 1);
    EXPECT_EQ(7, w.sum());
    EXPECT_EQ(4, w.elltwosquared());
    EXPECT_EQ(2, w.elltwo());
    EXPECT_EQ(hash, w.hash());
    // Test "merge" constructor
    Wheelwright::impl::Vector<int,Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> z(vals.begin(), vals.end(), 2, 3, 1, 1);
    Wheelwright::impl::Vector<int,Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> x(v, z);
    EXPECT_EQ(w.sum(), x.sum());
    EXPECT_EQ(w.elltwosquared(), x.elltwosquared());
    EXPECT_EQ(w.elltwo(), x.elltwo());
    Wheelwright::impl::Vector<int,Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> x2(z, v);
    EXPECT_EQ(x.sum(), x2.sum());
    EXPECT_EQ(x.elltwosquared(), x2.elltwosquared());
    EXPECT_EQ(x.elltwo(), x2.elltwo());
    EXPECT_EQ(x.hash(), x2.hash());
    // Test operator<
    EXPECT_FALSE(v < v);
    // Test operator==
    EXPECT_TRUE(v == v);
    EXPECT_FALSE(v == w);
    // Test that the computation of the hash satisfies basic requirements for
    // vectors that are prefixes/suffixes of others.
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> a(
                vals.begin(), vals.end(), 1, 0, 1);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> b(
                vals.begin(), vals.end(), 1, 1, 1);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA|Wheelwright::Features::DISCRETE> c(
                vals.begin(), vals.end(), 1, 0, 2);
    EXPECT_FALSE(a.hash() == b.hash());
    EXPECT_FALSE(a.hash() == c.hash());
    EXPECT_FALSE(b.hash() == c.hash());
}

TEST(VectorTest, VectorVANILLAONEDRAWWorks) {
    // Test constructor with values (and sum(), elltwosquared(), hash(), and
    // elltwo())
    std::size_t hash {0};
    boost::hash_combine(hash, 0u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 1.0);
    std::vector<Wheelwright::impl::Rademacher> pos(1,
            Wheelwright::impl::Rademacher::POS);
    std::vector<Wheelwright::impl::Rademacher> neg(1,
            Wheelwright::impl::Rademacher::NEG);
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> v(1, 0,
                1, neg);
    EXPECT_EQ(1, v.sum());
    EXPECT_EQ(1, v.elltwosquared());
    EXPECT_EQ(1, v.elltwo());
    EXPECT_EQ(-1, v.onedrawsum());
    EXPECT_EQ(hash, v.hash());
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> h(1, 0,
                1, pos, 2);
    EXPECT_EQ(1, h.sum());
    EXPECT_EQ(4, h.elltwosquared());
    EXPECT_EQ(2, h.elltwo());
    EXPECT_EQ(1, h.onedrawsum());
    EXPECT_EQ(v.hash(), h.hash());
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> h2(1, 0,
                1, pos, std::numeric_limits<double>::max());
    EXPECT_EQ(1, h2.sum());
    EXPECT_EQ(1, h2.elltwosquared());
    EXPECT_EQ(1, h2.elltwo());
    EXPECT_EQ(1, h2.onedrawsum());
    EXPECT_EQ(v.hash(), h.hash());
    // Test insert()
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2u);
    boost::hash_combine(hash, 2.0);
    std::vector<Wheelwright::impl::Rademacher> pos2(2, Wheelwright::impl::Rademacher::POS);
    v.insert(2, 1, 2, pos2);
    EXPECT_EQ(5, v.sum());
    EXPECT_EQ(9, v.elltwosquared());
    EXPECT_EQ(3, v.elltwo());
    EXPECT_EQ(3, v.onedrawsum());
    EXPECT_EQ(hash, v.hash());
    h.insert(2, 1, 2, pos2, 1);
    EXPECT_EQ(5, h.sum());
    EXPECT_EQ(6, h.elltwosquared());
    EXPECT_EQ(std::sqrt(6), h.elltwo());
    EXPECT_EQ(5, h.onedrawsum());
    h2.insert(2, 1, 2, pos2, std::numeric_limits<double>::max());
    EXPECT_EQ(5, h2.sum());
    EXPECT_EQ(9, h2.elltwosquared());
    EXPECT_EQ(3, h2.elltwo());
    EXPECT_EQ(5, h2.onedrawsum());
    // Test constructor with vector and new components
    boost::hash_combine(hash, 3u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2.0);
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> w(v, 2,
                3, 1, neg);
    EXPECT_EQ(7, w.sum());
    EXPECT_EQ(13, w.elltwosquared());
    EXPECT_EQ(std::sqrt(13), w.elltwo());
    EXPECT_EQ(hash, w.hash());
    EXPECT_EQ(1, w.onedrawsum());
    // Test "merge" constructor
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> z(2, 3,
                1, neg);
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> x(v, z);
    EXPECT_EQ(w.sum(), x.sum());
    EXPECT_EQ(w.elltwosquared(), x.elltwosquared());
    EXPECT_EQ(w.elltwo(), x.elltwo());
    EXPECT_EQ(w.onedrawsum(), x.onedrawsum());
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> x2(z, v);
    EXPECT_EQ(x.sum(), x2.sum());
    EXPECT_EQ(x.elltwosquared(), x2.elltwosquared());
    EXPECT_EQ(x.elltwo(), x2.elltwo());
    EXPECT_EQ(x.hash(), x2.hash());
    EXPECT_EQ(x.onedrawsum(), x2.onedrawsum());
    // Test operator<
    EXPECT_FALSE(v < v);
    // Test operator==
    EXPECT_TRUE(v == v);
    EXPECT_FALSE(v == w);
    // Test that the computation of the hash satisfies basic requirements for
    // vectors that are prefixes/suffixes of others.
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> a(1, 0,
                1, pos);
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> b(1, 1,
                1, pos);
    Wheelwright::impl::Vector<double,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::VANILLA> c(1, 0,
                2, pos2);
    EXPECT_FALSE(a.hash() == b.hash());
    EXPECT_FALSE(a.hash() == c.hash());
    EXPECT_FALSE(b.hash() == c.hash());
}

TEST(VectorTest, VectorONEDRAWDISCRETEWorks) {
    // Test constructor with values (and sum() and hash())
    std::size_t hash {0};
    std::set<int> vals {0, 1, 2, 3};
    boost::hash_combine(hash, 0u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 1);
    std::vector<Wheelwright::impl::Rademacher> pos(1,
            Wheelwright::impl::Rademacher::POS);
    std::vector<Wheelwright::impl::Rademacher> neg(1,
            Wheelwright::impl::Rademacher::NEG);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE>
            v(vals.begin(), vals.end(), 1, 0, 1, neg);
    EXPECT_EQ(1, v.sum());
    EXPECT_EQ(-1, v.onedrawsum());
    EXPECT_EQ(hash, v.hash());
    // Test domain_begin() and domain_end()
    std::size_t i {0};
    bool found_one {false};
    bool found_three {false};
    bool found_two {false};
    for (auto it {v.domain_begin()}; it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test insert()
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2u);
    boost::hash_combine(hash, 2);
    std::vector<Wheelwright::impl::Rademacher> pos2(2, Wheelwright::impl::Rademacher::POS);
    v.insert(2, 1, 2, pos2);
    EXPECT_EQ(5, v.sum());
    EXPECT_EQ(3, v.onedrawsum());
    EXPECT_EQ(hash, v.hash());
    // Test counts_begin() and counts_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.counts_begin(); it != v.counts_end(); ++it) {
        switch(*it) {
            case 0:
                EXPECT_FALSE(found_three);
                found_three = true;
                break;
            case 1:
                EXPECT_FALSE(found_one);
                found_one = true;
                break;
            case 2:
                EXPECT_FALSE(found_two);
                found_two = true;
                break;
            default:
                EXPECT_FALSE(true);
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test domain_begin() and domain_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.domain_begin(); it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test constructor with vector and new components
    boost::hash_combine(hash, 3u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE> w(v, 2,
                3, 1, neg);
    EXPECT_EQ(7, w.sum());
    EXPECT_EQ(1, w.onedrawsum());
    EXPECT_EQ(hash, w.hash());
    // Test "merge" constructor
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE>
            z(vals.begin(), vals.end(), 2, 3, 1, neg);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE> x(v, z);
    EXPECT_EQ(w.sum(), x.sum());
    EXPECT_EQ(w.onedrawsum(), x.onedrawsum());
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE> x2(z, v);
    EXPECT_EQ(x.sum(), x2.sum());
    EXPECT_EQ(x.onedrawsum(), x2.onedrawsum());
    EXPECT_EQ(x.hash(), x2.hash());
    // Test operator<
    EXPECT_FALSE(v < v);
    // Test operator==
    EXPECT_TRUE(v == v);
    EXPECT_FALSE(v == w);
    // Test that the computation of the hash satisfies basic requirements for
    // vectors that are prefixes/suffixes of others.
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE>
            a(vals.begin(), vals.end(), 1, 0, 1, pos);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE>
            b(vals.begin(), vals.end(), 1, 1, 1, pos);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE>
            c(vals.begin(), vals.end(), 1, 0, 2, pos2);
    EXPECT_FALSE(a.hash() == b.hash());
    EXPECT_FALSE(a.hash() == c.hash());
    EXPECT_FALSE(b.hash() == c.hash());
}

TEST(VectorTest, VectorONEDRAWDISCRETEVANILLAWorks) {
    // Test constructor with values (and sum(), elltwosquared(), hash(), and
    // elltwo())
    std::size_t hash {0};
    std::set<int> vals {0, 1, 2, 3};
    boost::hash_combine(hash, 0u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 1);
    std::vector<Wheelwright::impl::Rademacher> pos(1,
            Wheelwright::impl::Rademacher::POS);
    std::vector<Wheelwright::impl::Rademacher> neg(1,
            Wheelwright::impl::Rademacher::NEG);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA|Wheelwright::Features::ONEDRAW|
            Wheelwright::Features::DISCRETE> v(
                    vals.begin(), vals.end(), 1, 0, 1, neg, 0.5);
    EXPECT_EQ(1, v.sum());
    EXPECT_EQ(-1, v.onedrawsum());
    EXPECT_DOUBLE_EQ(0.25, v.elltwosquared());
    EXPECT_DOUBLE_EQ(0.5, v.elltwo());
    EXPECT_EQ(hash, v.hash());
    // Test domain_begin() and domain_end()
    std::size_t i {0};
    bool found_one {false};
    bool found_three {false};
    bool found_two {false};
    for (auto it {v.domain_begin()}; it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test insert()
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2u);
    boost::hash_combine(hash, 2);
    std::vector<Wheelwright::impl::Rademacher> pos2(2,
            Wheelwright::impl::Rademacher::POS);
    v.insert(2, 1, 2, pos2, 1);
    EXPECT_EQ(5, v.sum());
    EXPECT_EQ(3, v.onedrawsum());
    EXPECT_DOUBLE_EQ(2.25, v.elltwosquared());
    EXPECT_DOUBLE_EQ(std::sqrt(2.25), v.elltwo());
    EXPECT_EQ(hash, v.hash());
    // Test counts_begin() and counts_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.counts_begin(); it != v.counts_end(); ++it) {
        switch(*it) {
            case 0:
                EXPECT_FALSE(found_three);
                found_three = true;
                break;
            case 1:
                EXPECT_FALSE(found_one);
                found_one = true;
                break;
            case 2:
                EXPECT_FALSE(found_two);
                found_two = true;
                break;
            default:
                EXPECT_FALSE(true);
        }
        ++i;
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test domain_begin() and domain_end()
    i = 0;
    found_one = false;
    found_three = false;
    found_two = false;
    for (auto it = v.domain_begin(); it != v.domain_end(); ++it) {
        switch(*it) {
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
        }
        ++i;
    }
    EXPECT_TRUE(found_three);
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_two);
    EXPECT_EQ(3u, i);
    // Test constructor with vector and new components
    boost::hash_combine(hash, 3u);
    boost::hash_combine(hash, 1u);
    boost::hash_combine(hash, 2);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA|Wheelwright::Features::ONEDRAW|
            Wheelwright::Features::DISCRETE> w(v, 2, 3, 1, neg);
    EXPECT_EQ(7, w.sum());
    EXPECT_EQ(1, w.onedrawsum());
    EXPECT_DOUBLE_EQ(6.25, w.elltwosquared());
    EXPECT_DOUBLE_EQ(std::sqrt(6.25), w.elltwo());
    EXPECT_EQ(hash, w.hash());
    // Test "merge" constructor
    Wheelwright::impl::Vector<int,Wheelwright::Features::VANILLA|
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE>
        z(vals.begin(), vals.end(), 2, 3, 1, neg);
    Wheelwright::impl::Vector<int,Wheelwright::Features::VANILLA|
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE> x(v, z);
    EXPECT_EQ(w.sum(), x.sum());
    EXPECT_EQ(w.elltwosquared(), x.elltwosquared());
    EXPECT_EQ(w.elltwo(), x.elltwo());
    EXPECT_EQ(w.onedrawsum(), x.onedrawsum());
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA| Wheelwright::Features::ONEDRAW|
        Wheelwright::Features::DISCRETE> x2(z, v);
    EXPECT_EQ(x.sum(), x2.sum());
    EXPECT_EQ(x.elltwosquared(), x2.elltwosquared());
    EXPECT_EQ(x.elltwo(), x2.elltwo());
    EXPECT_EQ(x.onedrawsum(), x2.onedrawsum());
    EXPECT_EQ(x.hash(), x2.hash());
    // Test operator<
    EXPECT_FALSE(v < v);
    // Test operator==
    EXPECT_TRUE(v == v);
    EXPECT_FALSE(v == w);
    // Test that the computation of the hash satisfies basic requirements for
    // vectors that are prefixes/suffixes of others.
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA| Wheelwright::Features::ONEDRAW|
            Wheelwright::Features::DISCRETE> a(vals.begin(), vals.end(), 1, 0,
                    1, pos);
    Wheelwright::impl::Vector<int,
        Wheelwright::Features::VANILLA|Wheelwright::Features::ONEDRAW|
            Wheelwright::Features::DISCRETE> b(vals.begin(), vals.end(), 1, 1,
                    1, pos);
    Wheelwright::impl::Vector<int,Wheelwright::Features::VANILLA|
        Wheelwright::Features::ONEDRAW|Wheelwright::Features::DISCRETE>
            c(vals.begin(), vals.end(), 1, 0, 2, pos2);
    EXPECT_FALSE(a.hash() == b.hash());
    EXPECT_FALSE(a.hash() == c.hash());
    EXPECT_FALSE(b.hash() == c.hash());
}
