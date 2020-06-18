/** @file PtrLess.h
 * Like std::less<T>, but compare the objects pointed to by two shared_ptr's.
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

#ifndef MATTEOUTILS_PTRLESS_H_
#define MATTEOUTILS_PTRLESS_H_

#include <functional>
#include <memory>

namespace MatteoUtils {

/** Like std::less<T>, but compares the objects pointed to by two shared_ptr's.
 *
 * @tparam T     the type of the objects stored in the pointers.
 * @tparam Less  the functor to use to compare the objects. Default:
 *               std::less<T>.
 */
template <class T, class Less = std::less<T>> class PtrLess {
    public:
        /** Performs the comparison by calling @a Less with  @a a and @a b.
         *
         * @param a,b  the std::shared_ptr's pointing to the objects to compare.
         *
         * @returns  Less(a)(*a, *b) if @a a and @a b are not nullptr. true if
         *           @a a is nullptr. false if @a b is nullptr and @a is not
         *           nullptr.
         */
        bool operator() (const std::shared_ptr<T> &a,
                const std::shared_ptr<T> &b) const {
            if (a == nullptr)
                return b != nullptr;
            else if (b == nullptr)
                return false;
            else
                return Less()(*a, *b);
        }
}; // class PtrLess

} // namespace MatteoUtils

#endif // MATTEOUTILS_PTRLESS_H_
