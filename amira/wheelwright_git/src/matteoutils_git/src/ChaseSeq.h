/** @file ChaseSeq.h
 * Generate combinations of objects in Chase's sequence order.
 *
 * @author Matteo Riondato
 * @date 2017 06 07
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

#ifndef MATTEOUTILS_CHASESEQ_H_
#define MATTEOUTILS_CHASESEQ_H_

#include <cstddef>
#include <iterator>
#ifdef MATTEOUTILS_INVARIANTS
#include <stdexcept>
#endif

#include "ChaseSeqIter.h"
#include "MissingMath.h"

namespace MatteoUtils {

/** Generate combinations of objects in Chase's sequence order.
 *
 * This class is mostly an interface to the ChaseSeqIter class.
 *
 * @tparam T  the type of the elements whose combinations to generate.
 */
template<class T> class ChaseSeq {
    public:
        /** The type returned by begin() and end(). */
        using iterator = impl::ChaseSeqIterator<T>;

        /** The type returned by cbegin() and cend(). */
        using const_iterator = impl::ChaseSeqConstIterator<T>;

        /** Constructor.
         *
         * @tparam  Iterator  the type of the elements.
         *
         * @param [begin,end)  the sequence of the elements.
         * @param maxl         the maximum length of the combinations.
         * @param minl         the minimum length of the combinations.
         */
        template<class Iterator> ChaseSeq(const Iterator &begin,
                const Iterator &end, const std::size_t maxl,
                const std::size_t minl = 1) : maxl {maxl}, minl {minl},
                elems {begin, end},
                sz {sum_of_binoms(elems.size(), maxl, minl)} {
            const std::set<T> s(elems.cbegin(), elems.cend());
#ifdef MATTEOUTILS_INVARIANTS
            if (maxl < minl)
                throw std::runtime_error("MatteoUtils::ChaseSeq::ChaseSeq: the "
                        "max requested length cannot be smaller than the min "
                        "one.");
            if (s.size() != elems.size())
                throw std::runtime_error("MatteoUtils::ChaseSeq::ChaseSeq: the "
                        "container from 'begin' to 'end' cannot contain "
                        "repeated elements.");
#endif
        }

        /** Return an iterator to the first combination. */
        iterator begin() const {
            return iterator(elems, maxl, minl);
        }

        /** Return an iterator past the last combination. */
        iterator end() const {
            return iterator(elems, maxl, minl, impl::ChaseSeqIterState::END);
        }

        /** Return a const iterator to the first combination. */
        const_iterator cbegin() const {
            return const_iterator(elems, maxl, minl);
        }

        /** Return a const iterator past the last combination. */
        const_iterator cend() const {
            return const_iterator(elems, maxl, minl,
                    impl::ChaseSeqIterState::END);
        }

        /** Return the number of combinations in the sequence. */
        std::size_t size() const { return sz; }

    private:
        const std::size_t maxl;
        const std::size_t minl {1};
        const std::vector<T> elems;
        const std::size_t sz;
}; // class ChaseSeq

} // namespace MatteoUtils

#endif // MATTEOUTILS_CHASESEQ_H_
