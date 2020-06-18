/** @file ChaseSeqIter.h
 * Generate combinations of elements in Chase's sequence order.
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

#ifndef MATTEOUTILS_CHASESEQITER_H_
#define MATTEOUTILS_CHASESEQITER_H_

#include <cassert>
#include <set>
#ifdef MATTEOUTILS_INVARIANTS
#include <stdexcept>
#endif
#include <utility>
#include <vector>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include "MissingMath.h"

namespace MatteoUtils {

template<class T> class ChaseSeq; /* forward declaration */

} // namespace MatteoUtils

namespace MatteoUtils::impl {

/** Initial state of the iterator */
enum class ChaseSeqIterState : char {
    /** begin */
    BEGIN = 1,
    /** end */
    END   = 2
};

/** Generate combinations of elements in a set in Chase's sequence order.
 *
 * This class implements a straightforward extension to Algorithm C (Chase's
 * sequence) from D. Knuth's TAOCP, Vol. 4A, Sect. 7.2.1.3, page 367. The
 * extension consists in allowing to generate all combinations of elements from
 * a set, from a minimum combination length to a maximum combination length
 * (extremes included).
 *
 * It is intended to be used only through the ChaseSeq class.
 *
 * @tparam T      the type of the elements.
 * @tparam Value  the type of the 'pointed-to' objects, i.e., (const)
 *                std::pair<const std::set<T>::const_iterator,const
 *                std::set<T>::const_iterator>.
 */
template<class T, class Value> class ChaseSeqIter :
        public boost::iterator_facade<ChaseSeqIter<T,Value>, Value,
        boost::single_pass_traversal_tag, Value> {
    friend class boost::iterator_core_access;
    friend class ::MatteoUtils::ChaseSeq<T>;

    template <class OtherValue> bool equal(
            ChaseSeqIter<T,OtherValue> const& o) const {
        return (minl == o.minl && maxl == o.maxl && elems == o.elems
                && progr == o.progr);
    }

    Value dereference() const {
        return std::make_pair(comb.cbegin(), comb.cend());
    }

    void increment() { c3(); ++progr; }

    struct enabler {}; // just for boost::enable_if

    struct in_out {
        std::size_t in;
        std::size_t out;
    };

    /* The names of the following private member functions correspond to the
     * names of the steps in the pseudocode in TAOCP. The only difference is
     * step c2, that here corresponds to dereference and increment.
     */
    void c3();
    void c4(in_out &);
    void c5(in_out &);
    void c6(in_out &);
    void c7(in_out &);
    void reset();

    bool end {false};
    bool endt {false};
    std::size_t j {0};
    std::size_t progr {0};
    const std::size_t minl;
    const std::size_t maxl;
    const std::vector<T> &elems;
    std::size_t r;
    std::size_t l;
    std::vector<bool> a;
    std::vector<bool> w;
    std::set<T> comb;

    /** Constructor.
     *
     * @param elems  reference to a std::vector<T> of elements. The vector
     *               should actually represent a set, i.e., the elements must
     *               be distinct, but a vector allows for random access, which
     *               we need.
     * @param maxl   the maximum combination length.
     * @param minl   the minimum combination length.
     * @param state  the initial state of the iterator.
     */
    ChaseSeqIter(const std::vector<T> &elems, const std::size_t maxl,
            const std::size_t minl = 1,
            const ChaseSeqIterState state = ChaseSeqIterState::BEGIN);

    /** Converting constructor.
     *
     * @tparam OtherValue  the type for the second template parameter of
     *                     the ChaseSeqIter<l,OtherValue> object to `copy''.
     *
     * @param o  the ChaseSeqIter<l,OtherValue> object to ``copy''.
     */
    template <class OtherValue>
            ChaseSeqIter(ChaseSeqIter<T, OtherValue> const& o,
                typename boost::enable_if<
                boost::is_convertible<OtherValue*,Value*>,enabler>::type =
                    enabler()) :
                end {o.end}, endt {o.endt}, j {o.j}, progr {o.progr},
                minl {o.minl}, maxl {o.maxl}, elems {o.elemems_v},
                r {o.r}, l {o.l}, a {o.a}, w {o.w}, comb {o.comb} {}

}; // class ChaseSeqIter


template<class T, class Value> ChaseSeqIter<T,Value>::ChaseSeqIter(
        const std::vector<T> &elems, const std::size_t maxlen,
        const std::size_t minlen, ChaseSeqIterState state) :
        progr {(state == ChaseSeqIterState::BEGIN) ? 0 :
            sum_of_binoms(elems.size(), maxlen, minlen)},
        minl {minlen}, maxl {maxlen}, elems {elems},
        r {(elems.size() > minlen)? elems.size() - minlen : minlen},
        l {minlen}, a(elems.size(), false), w (elems.size() + 1, true) {
#ifdef MATTEOUTILS_INVARIANTS
    if (minlen > maxlen)
        throw std::invalid_argument("MatteoUtils::ChaseSeqIter::ChaseSeqIter: "
                "minlen must not be larger than maxlen.");
    if (maxlen > elems.size())
        throw std::invalid_argument("MatteoUtils::ChaseSeqIter::ChaseSeqIter: "
                "maxlen must not be larger than the number of elements.");
#endif
    a.shrink_to_fit();
    w.shrink_to_fit();
    for (j = elems.size() - l; j < elems.size(); ++j) {
        a[j] = true;
        comb.insert(elems[j]);
    }
}


template<class T, class Value> void ChaseSeqIter<T,Value>::c3() {
    assert(! endt && ! end);
    for (j = r; ! w[j]; ++j)
        w[j] = true;
    if (j == elems.size()) {
        endt = true;
        if (l < maxl) {
            ++l;
            reset();
        } else
            end = true;
    } else {
        assert(w[j]);
        w[j] = false;
        in_out chg;
        if (a[j]) {
            if (j % 2 == 1)
                c4(chg);
            else
                c5(chg);
        } else {
            if (j % 2 == 1)
                c7(chg);
            else
                c6(chg);
        }
        assert(! a[chg.in] && a[chg.out]);
        a[chg.in] = true;
        a[chg.out] = false;
        assert(comb.count(elems[chg.out]) == 1);
        assert(comb.count(elems[chg.in]) == 0);
        comb.erase(elems[chg.out]);
#ifdef NDEBUG
        comb.insert(elems[chg.in]);
#else
        assert(comb.insert(elems[chg.in]).second);
#endif
    }
}


template<class T, class Value> void ChaseSeqIter<T,Value>::c4(in_out &p) {
    assert(! endt && ! end);
    assert(j >= 1);
    p.in = j - 1;
    p.out = j;
    if (r == j && j > 1)
        r = j - 1;
    else if (r == j - 1)
        r = j;
}


template<class T, class Value> void ChaseSeqIter<T,Value>::c5(in_out &p) {
    assert(! endt && ! end);
    assert(j >= 2);
    if (a[j - 2])
        c4(p);
    else {
        p.in = j - 2;
        p.out = j;
        if (r == j)
            r = std::max(j - 2, (std::size_t) 1);
        else if (r == j - 2)
            r = j - 1;
    }
}


template<class T, class Value> void ChaseSeqIter<T,Value>::c6(in_out &p) {
    assert(! endt && ! end);
    assert(j > 0);
    p.in = j;
    p.out = j - 1;
    if (r == j && j > 1)
        r = j - 1;
    else if (r == j - 1)
        r = j;
}


template<class T, class Value> void ChaseSeqIter<T,Value>::c7(in_out &p) {
    assert(! endt && ! end);
    assert(j > 0);
    if (a[j - 1])
        c6(p);
    else {
        assert(j >= 2);
        p.in = j;
        p.out = j - 2;
        if (r == j - 2)
            r = j;
        else if (r == j - 1)
            r = j - 2;
    }
}


template<class T, class Value> void ChaseSeqIter<T,Value>::reset() {
    endt = false;
    end = false;
    comb.clear();
    a.assign(elems.size(), false);
    for (j = elems.size() - l; j < elems.size(); ++j) {
        a[j] = true;
        comb.insert(elems[j]);
    }
    w.assign(elems.size() + 1, true);
    r = (elems.size() > l) ? elems.size() - l : l;
}

/** Useful alias. */
template<class T> using setcit = typename std::set<T>::const_iterator;

/** Useful alias. */
template<class T> using ChaseSeqIterator =
ChaseSeqIter<T,std::pair<setcit<T>,setcit<T>>>;

/** Useful alias. */
template<class T> using ChaseSeqConstIterator =
ChaseSeqIter<T, const std::pair<setcit<T>,setcit<T>>>;

} // namespace MatteoUtils::impl

#endif // MATTEOUTILS_CHASESEQITER_H_
