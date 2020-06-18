/** @file PtrHash.h
 * Like std::hash<T>, but computes the hash for an object pointed to by a
 * std::shared_ptr.
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

#ifndef MATTEOUTILS_PTRHASH_H_
#define MATTEOUTILS_PTRHASH_H_

#include <functional>
#include <memory>
#include <stdexcept>

namespace MatteoUtils {

/** Like std::hash<T>, but computes the hash for an object pointed to by a
 * std::shared_ptr.
 *
 * @tparam T     the type of the objects stored in the std::shaared_ptr.
 * @tparam Hash  the function class to use to compute the hash. Default:
 *               std::hash<T>.
 */
template <class T, class Hash = std::hash<T>> class PtrHash {
    public:
        /** Computes the hash by calling a @a Hash functor on @a a.
         *
         * @param a  the object whose hash to compute.
         */
        std::size_t operator() (const std::shared_ptr<T> &a) const {
#ifdef MATTEOUTILS_INVARIANTS
            if (a == nullptr)
                // TODO:  MatteoUtils::PtrHash : what do do in case of nullptr? It
                // seems that in C++1z there exists a specialization
                // std::hash<std::nullptr_t>, but how to use it?
                throw std::invalid_argument("MatteoUtils::PtrHash::operator(): "
                        "the passed pointer equals nullptr.");
#endif
            return Hash{}(*a);
        }
}; // class PtrHash

} // namespace MatteoUtils

#endif // MATTEOUTILS_PTRHASH_H_
