/** @file Iterators.h
 * Template metaprogramming functions to check iterators
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

#ifndef MATTEOUTILS_ITERATORS_H
#define MATTEOUTILS_ITERATORS_H

#include <iterator>
#include <type_traits>

namespace MatteoUtils {

/** Check whether the type @a T is an iterator.
 *
 * Implementation from https://stackoverflow.com/a/12032923 .
 *
 * @tparam T  the type to check.
 */
template<typename T, typename = void> struct is_iterator {
    static constexpr bool value = false;
};

/** Check whether the type @a T is an iterator.
 *
 * Implementation from https://stackoverflow.com/a/12032923 .
 *
 * @tparam T  the type to check.
 */
template<typename T>
struct is_iterator<T, typename std::enable_if<!std::is_same<typename
std::iterator_traits<T>::value_type, void>::value>::type> {
    static constexpr bool value = true;
};

/** Check whether the type @a T is a random access iterator.
 *
 * @tparam T  the type to check.
 */
template<typename T, typename = void> struct is_random_access_iterator {
    static constexpr bool value = false;
};

/** Check whether the type @a T is a random access iterator.
 *
 * @tparam T  the type to check.
 */
template<typename T>
struct is_random_access_iterator<T,
    typename std::enable_if<std::is_same<typename std::iterator_traits<T>::iterator_category, std::random_access_iterator_tag>::value>::type> {
    static constexpr bool value = true;
};

} // namespace MatteoUtils
#endif // MATTEOUTILS_ITERATORS_H
