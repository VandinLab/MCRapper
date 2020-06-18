/** @file RandomAccessIteratorIncludes.h
 * A specialization of std::include for random access iterators.
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

#ifndef MATTEOUTILS_RANDOMACCESSITERATORINCLUDES_H_
#define MATTEOUTILS_RANDOMACCESSITERATORINCLUDES_H_

#include <algorithm>
#include <iterator>

namespace MatteoUtils {

/** Checks whether a sorted sequence is a subsequence of another.
 *
 *  @tparam T1  the type of the iterator for the first sequence.
 *  @tparam T2  the type of the iterator for the second sequence.
 *
 *  @param first1  iterator to the first element of the first sequence.
 *  @param last1   iterator past the last element of the first sequence.
 *  @param first2  iterator to the first element of the second sequence.
 *  @param last2   iterator past the last element of the second sequence.
 */
template<typename T1,typename T2> bool includes(T1 first1, T1 last1,
        T2 first2, T2 last2, std::enable_if_t<
         std::is_same_v<typename std::iterator_traits<T1>::iterator_category,
         std::random_access_iterator_tag> && std::is_same_v<typename
         std::iterator_traits<T2>::iterator_category,
         std::random_access_iterator_tag>,bool> = false) {
    while (std::distance(first1, last1) >= std::distance(first2, last2) &&
            first2 < last2) {
        if (auto b {std::lower_bound(first1, last1, *first2)};
                b != last1 && *b == *first2)
            first1 = b + 1;
        else
            break;
        ++first2;
    }
    return first2 == last2;
}

} // namespace MatteoUtils

#endif // MATTEOUTILS_RANDOMACCESSITERATORINCLUDES_H_
