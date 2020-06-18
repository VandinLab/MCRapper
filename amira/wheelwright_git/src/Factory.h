/** @file Factory.h
 * Factory for Projection vectors.
 *
 * @author Matteo Riondato
 * @date 2017 09 08
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

#ifndef WHEELWRIGHT_FACTORY_H
#define WHEELWRIGHT_FACTORY_H

#include <cstddef>
#include <memory>
#include <set>

#include "Features.h"

namespace Wheelwright {
    enum Rademacher : signed char;
} // namespace Wheelwright

namespace Wheelwright::impl {

template<typename T, int F, typename Less, typename J> class Vector;

/** Factory for Projection vectors, including the discrete variant.
 *
 * The function members of this class mimic the ones with the same arguments
 * from rojVec and DiscrVector.
 *
 * @tparam T the Vector type
 */
template<typename T> class Factory {
    public:
        std::shared_ptr<T> make(const T &e) {
            return std::make_shared<T>(e);
        }

        template<typename U = T> std::shared_ptr<T> make(const T &e,
                const typename T::value_type v, const std::size_t i,
                const std::size_t c, std::enable_if_t<U::features ==
                    Wheelwright::Features::DISCRETE> = 0) {
            return std::make_shared<T>(e, v, i, c);
        }

        template<typename U = T> std::shared_ptr<U> make(const U &e,
                const typename U::value_type v, const std::size_t i,
                const std::size_t c,
                const double lv = std::numeric_limits<double>::max()//,
                //std::enable_if_t<U::features ==
                //    Wheelwright::Features::VANILLA> = 0
                ) {
            return std::make_shared<U>(e, v, i, c, lv);
        }

        std::shared_ptr<T> make(const typename T::value_type v,
                const std::size_t i, const std::size_t c) {
            return std::make_shared<T>(v, i, c);
        }

        template<typename... S> std::shared_ptr<T> make(
                const typename T::value_type v, const std::size_t i,
                const std::size_t c, S... rest) {
            return std::make_shared<T>(v, i, c, rest...);
        }

        std::shared_ptr<T> make(const T &a, const T &b) {
            return std::make_shared<T>(a, b);
        }
}; // class Factory


template<typename T> class DiscreteFactoryBase {

    public:
        template<typename Iterator> DiscreteFactoryBase(Iterator b,
                Iterator e, typename std::enable_if_t<!
                std::is_arithmetic<Iterator>::value>* = 0) : vals(b, e) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (vals.size() == 0)
                throw std::runtime_error("Wheelwright::Factory: "
                        "vals cannot be empty.");
#endif
        }

        template<typename... J> DiscreteFactoryBase(J ...args) : vals {args...} {
#ifdef WHEELWRIGHT_INVARIANTS
            if (vals.size() == 0)
                throw std::runtime_error("Wheelwright::Factory: "
                        "vals cannot be empty.");
#endif
        }

    protected:
        template<int F, typename Less> std::shared_ptr<Vector<T,F,Less,T>> make(
                const Vector<T,F,Less,T> &e) {
            return std::make_shared<Vector<T,F,Less,T>>(e);
        }

        template<int F, typename Less> std::shared_ptr<Vector<T,F,Less,T>> make(
                const Vector<T, F, Less, T> &a,
                const Vector<T, F, Less, T> &b) {
            return std::make_shared<Vector<T,F,Less,T>>(a, b);
        }

        std::set<T> vals;
};

/** Specialization of the factory for discrete projection vectors. */
template<typename T, int F, typename Less> class
Factory<Vector<T,F,Less,
    std::enable_if_t<F == Features::DISCRETE, T>>> : DiscreteFactoryBase<T> {

    using parent = DiscreteFactoryBase<T>;
    using parent::vals;

    public:
        using parent::parent;
        using parent::make;

        std::shared_ptr<Vector<T,F,Less,T>> make(const Vector<T,F,Less,T> &e,
                const T v, const std::size_t i, const std::size_t c) {
            return std::make_shared<Vector<T,F,Less,T>>(e, v, i, c);
        }

        std::shared_ptr<Vector<T,F,Less,T>> make(
                const T v, const std::size_t i, const std::size_t c) {
            return std::make_shared<Vector<T,F,Less,T>>(vals.cbegin(),
                    vals.cend(), v, i, c);
        }
};

template<typename T, int F, typename Less> class
Factory<Vector<T,F,Less,
    std::enable_if_t<F == (Features::DISCRETE |
            Features::VANILLA), T>>> : DiscreteFactoryBase<T> {
    using parent = DiscreteFactoryBase<T>;
    using parent::vals;

    public:
        using parent::parent;
        using parent::make;

        std::shared_ptr<Vector<T,F,Less,T>> make(
                const Vector<T,F,Less,T> &e, const T v,
                const std::size_t i, const std::size_t c,
                const double lv = std::numeric_limits<double>::max()) {
            return std::make_shared<Vector<T,F,Less,T>>(e, v, i, c, lv);
        }

        std::shared_ptr<Vector<T,F,Less,T>> make(
                const T v, const std::size_t i, const std::size_t c,
                const double lv = std::numeric_limits<double>::max()) {
            return std::make_shared<Vector<T,F,Less,T>>(vals.cbegin(),
                        vals.cend(), v, i, c, lv);
        }
};

} // namespace Wheelwright::impl

#endif // WHEELWRIGHT_FACTORY_H
