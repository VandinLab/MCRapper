/** @file KeyValueGetters.h
 *
 * Unary functors to return the key or the value of an iterator of a std::map or
 * std::unordered_map.
 *
 * @author Matteo Riondato
 * @date 2017 09 10
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

#ifndef MATTEOUTILS_KEYVALUEGETTERS_H_
#define MATTEOUTILS_KEYVALUEGETTERS_H_

namespace MatteoUtils {

/** Returns the first element of a pair obtained from dereferencing an iterator
 * for a std::map or std::unordered_map.
 *
 * @tparam Iter  the type of the iterator.
 */
template<typename Iter> struct KeyGetter {
    /** Returns the first element of the pair @a p.  */
    const typename Iter::value_type::first_type& operator()
            (const typename Iter::value_type& p) const {
        return p.first;
    }
};

/** Returns the second element of a pair obtained from dereferencing an iterator
 * for a std::map or std::unordered_map.
 *
 * @tparam Iter  the type of the iterator.
 */
template<typename Iter> struct ValueGetter {
    /** Returns the second element of the pair @a p.  */
    const typename Iter::value_type::second_type& operator()
            (const typename Iter::value_type& p) const {
        return p.second;
    }
};

} // namespace MatteoUtils

#endif // MATTEOUTILS_KEYVALUEGETTERS_H_
