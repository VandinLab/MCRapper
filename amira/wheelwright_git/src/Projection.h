/** @file Projection.h
 * Projection set.
 *
 * @author Matteo Riondato
 * @date 2017 06 19
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

#ifndef WHEELWRIGHT_PROJECTION_H
#define WHEELWRIGHT_PROJECTION_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

//#include <boost/iterator/transform_iterator.hpp>

#include "matteoutils/KeyValueGetters.h"
#include "matteoutils/PtrEqualTo.h"
#include "matteoutils/PtrHash.h"
#include "matteoutils/PtrLess.h"

#include "Features.h"
#include "Projection_impl.h"
#include "Vector.h"

namespace Wheelwright {

/** Projection of a family of functions on to a sample.
 *
 * # Introduction
 *
 * Let \f$F\f$ be a family of functions from a domain \f$D\f$ to discrete finite
 * subset of reals. Assume that the functions in \f$F\f$ can be distinguished
 * from each other based on the value of a parameter \f$t\f$, belonging to some
 * domain \f$T\f$. I.e., \f$F=\{f_t, t\in T\}\f$.
 *
 * Let now \f$S=\{x_1,\dotsc,x_n\}\f$ be a sample from \f$D\f$, i.e., a bag of
 * \f$n\f$ elements from \f$D\f$. For each \f$t\in T\f$, let \f$v_t(S)\f$
 * be the \f$n\f$-dimensional vector \f$v_t(S)=\langle f_t(x_1), \dotsc,
 * f_t(x_n)\rangle\f$. The *projection* *of* \f$F\f$ *on* \f$S\f$ is the *set*
 * \f$P_F(S)\f$ of vectors \f$v_t(S)\f$, \f$t\in T\f$: \f$P_F(S)=\{v_t(S), t\in
 * T\}$. \f$P_F(S)\f$ is a *set* not a bag, so the number of vectors it
 * contains may be much smaller than the size of \f$T\f$.
 *
 * The average of the components of a vector \f$v\in P_F(S)\f$ gives the sample
 * average of the functions \f$f_t\in F\f$ such that \f$v_t=v\f$.
 *
 * The \f$\ell_2\f$-norm of the vectors in \f$P_F(S)\f$ is used to compute
 * probabilistic upper bounds to the empirical Rademacher averages and therefore
 * to the maximum deviation of the sample averages from their expectations.
 *
 * # Description of the class
 *
 * This class maintains \f$P_F(S)$ as the sample \f$S\f$ grows. It also
 * maintains a mapping from parameters in \f$T\f$ to vectors in \f$P_F(S)\f$.
 *
 * TODO: Wheelwright::Projection: update and complete documentation
 *
 * @tparam T        the type of the parameters mapped to vectors.
 * @tparam V        the type of the components of the vectors. Default: double.
 * @tparam Vector   the type of the vectors. See also the Vector<V> and
 *                  the DiscrVector<V> classes. Default: Vector<V>.
 * @tparam Hash     the type of a functor implementing a hash function of type
 *                  T. Default: std::hash<T>.
 * @tparam EqualTo  the type of a functor to test for equality of parameters of
 *                  type T. Default: std::equal_to<T>
 * @tparam LessV    the type of a functor implementing less-than for parameters
 *                  of type V. Default: std::less<V>.
 */
template <typename T,int F,typename V = double,typename Hash = std::hash<T>,
         typename EqualTo = std::equal_to<T>, typename LessV = std::less<V>>
         class Projection { };


template <typename T,typename V,typename Hash,typename EqualTo,typename LessV>
class Projection<T,Wheelwright::Features::VANILLA,V, Hash,EqualTo,LessV> :
public impl::ProjectionBase<T,Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV> {
    private:
        using Vector = impl::Vector<V,Wheelwright::Features::VANILLA,LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Easy alias for the ProjectionBase parent */
        using Parent = impl::ProjectionBase<T,Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>;
        using ThisClass = Projection<T,Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>;

        void add_zeroes(vec_p &v) override {
                v->insert(0, this->u_idx, this->u_cps, this->u_avg);
        }

        vec_p built_from_scratch(const V v) const override {
                return std::make_shared<Vector>(v, this->u_idx, this->u_cps,
                        v - this->u_avg);
        }

        vec_p extended(const Vector &vec, V v) const override {
            return std::make_shared<Vector>(vec, v, this->u_idx, this->u_cps,
                    v - this->u_avg);
        }

    public:
        explicit Projection(const std::size_t cr) : Parent {cr} {}

        /** Merge
         */
        ThisClass& merge(ThisClass& oth) {
            Parent::do_merge(oth);
            return *this;
        }

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.
}; // VANILLA


template <typename T,typename V,typename Hash,typename EqualTo,typename LessV>
class Projection<T,Wheelwright::Features::DISCRETE,V,Hash,EqualTo,LessV> :
public impl::ProjectionBase<T,Wheelwright::Features::DISCRETE,V,Hash,EqualTo,LessV>,
       public impl::DiscreteBase<V> {
    private:
        using Vector = impl::Vector<V,Wheelwright::Features::DISCRETE,LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Easy alias for the ProjectionBase parent */
        using Parent = impl::ProjectionBase<T,Wheelwright::Features::DISCRETE,V,Hash,EqualTo,LessV>;
        using ThisClass = Projection<T,Wheelwright::Features::DISCRETE,V,Hash,EqualTo,LessV>;

        void add_zeroes(vec_p &v) override {
                v->insert(0, this->u_idx, this->u_cps);
        }

        vec_p built_from_scratch(const V v) const override {
            return std::make_shared<Vector>(this->codomain.cbegin(),
                    this->codomain.cend(), v, this->u_idx, this->u_cps);
        }

        vec_p extended(const Vector &vec, V v) const override {
                return std::make_shared<Vector>(vec, v, this->u_idx,
                        this->u_cps);
        }

    public:
        /** Constructor.
         *
         * @tparam S  variadic template type for the types of the additional
         *            arguments to the constructor.
         *
         * @param cr    the number of references for the catch-all vector. I.e.,
         *              the size of the family of functions, if finite. An
         *              infinite function can be partially simulated by passing
         *              std::numeric_limits<T>::max().
         * @param rest  additional arguments forwarded to the constructore of a
         *              VectorFactory<Vector> object.
         */
        template<typename... S> explicit Projection(const std::size_t cr,
                S... rest) : Parent {cr}, impl::DiscreteBase<V> {rest...} {}

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        /** Merge
         */
        ThisClass& merge(ThisClass &oth) {
            impl::DiscreteBase<V>::do_merge(oth);
            Parent::do_merge(oth);
            return *this;
        }
}; // DISCRETE


template <typename T,typename V,typename Hash,typename EqualTo,typename LessV>
class Projection<T,Wheelwright::Features::DISCRETE |
Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV> :
public impl::ProjectionBase<T,Wheelwright::Features::DISCRETE |
Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>,
public impl::DiscreteBase<V> {
    private:
        using Vector = impl::Vector<V,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::VANILLA,LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Easy alias for the ProjectionBase parent */
        using Parent = impl::ProjectionBase<T,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>;
        using ThisClass = Projection<T,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>;

        void add_zeroes(vec_p &v) override {
                v->insert(0, this->u_idx, this->u_cps, this->u_avg);
        }

        vec_p built_from_scratch(const V v) const override {
                return std::make_shared<Vector>(this->codomain.cbegin(),
                        this->codomain.cend(), v, this->u_idx, this->u_cps,
                        v - this->u_avg);
        }

        vec_p extended(const Vector &vec, V v) const override {
                return std::make_shared<Vector>(vec, v, this->u_idx,
                        this->u_cps, v - this->u_avg);
        }

    public:
        /** Constructor.
         *
         * @tparam S  variadic template type for the types of the additional
         *            arguments to the constructor.
         *
         * @param cr    the number of references for the catch-all vector. I.e.,
         *              the size of the family of functions, if finite. An
         *              infinite function can be partially simulated by passing
         *              std::numeric_limits<T>::max().
         * @param rest  additional arguments forwarded to the constructore of a
         *              VectorFactory<Vector> object.
         */
        template<typename... S> explicit Projection(const std::size_t cr,
                S... rest) : Parent {cr}, impl::DiscreteBase<V> {rest...} {}

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        /** Merge
         */
        ThisClass& merge(ThisClass &oth) {
            impl::DiscreteBase<V>::do_merge(oth);
            Parent::do_merge(oth);
            return *this;
        }
}; // DISCRETE | VANILLA


template <typename T,typename V,typename Hash,typename EqualTo,typename LessV>
class Projection<T,Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV> :
public impl::ProjectionBase<T,Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV>,
       public impl::OneDrawBase {
    private:
        using Vector = impl::Vector<V,Wheelwright::Features::ONEDRAW,LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Easy alias for the ProjectionBase parent */
        using Parent = impl::ProjectionBase<T,Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV>;
        using ThisClass= Projection<T,Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV>;

        void add_zeroes(vec_p &vec) override {
            vec->insert(0, this->u_idx, this->u_cps, u_rade);
        }

        vec_p built_from_scratch(const V v) const override {
            return std::make_shared<Vector>(v, this->u_idx, this->u_cps,
                    u_rade);
        }

        vec_p extended(const Vector &vec, V v) const override {
            return std::make_shared<Vector>(vec, v, this->u_idx, this->u_cps,
                    u_rade);
        }

    public:
        explicit Projection(const std::size_t cr) :
            impl::ProjectionBase<
            T,Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV> {cr},
            impl::OneDrawBase {} {}

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        void begin_update(const std::size_t idx, const std::size_t cps) override {
            sample_rademacher(cps);
            Parent::begin_update(idx, cps);

        }

        void clear() override {
            Parent::clear();
            u_rade.clear();
        }

        void end_update() override {
            Parent::end_update();
            u_rade.clear();
        }

        /** Merge
         */
        ThisClass& merge(ThisClass& oth) {
            Parent::do_merge(oth);
            return *this;
        }
}; // ONEDRAW


template <typename T,typename V,typename Hash,typename EqualTo,typename LessV>
class Projection<T,Wheelwright::Features::DISCRETE |
Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV> : public impl::ProjectionBase<T,
    Wheelwright::Features::DISCRETE | Wheelwright::Features::ONEDRAW,V,Hash,
    EqualTo,LessV>, public impl::DiscreteBase<V>, public impl::OneDrawBase {
    private:
        using Vector = impl::Vector<V,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW,LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Easy alias for the ProjectionBase parent */
        using Parent = impl::ProjectionBase<T,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV>;
        using ThisClass = Projection<T,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW,V,Hash,EqualTo,LessV>;

        void add_zeroes(vec_p &vec) override {
            vec->insert(0, this->u_idx, this->u_cps, u_rade);
        }

        vec_p built_from_scratch(const V v) const override {
            return std::make_shared<Vector>(this->codomain.cbegin(),
                    this->codomain.cend(), v, this->u_idx, this->u_cps, u_rade);
        }

        vec_p extended(const Vector &vec, V v) const override {
            return std::make_shared<Vector>(vec, v, this->u_idx, this->u_cps,
                    u_rade);
        }

    public:
        /** Constructor.
         *
         * @tparam S  variadic template type for the types of the additional
         *            arguments to the constructor.
         *
         * @param cr    the number of references for the catch-all vector. I.e.,
         *              the size of the family of functions, if finite. An
         *              infinite function can be partially simulated by passing
         *              std::numeric_limits<T>::max().
         * @param rest  additional arguments forwarded to the constructore of a
         *              VectorFactory<Vector> object.
         */
        template<typename... S> explicit Projection(const std::size_t cr,
                S... rest) : Parent {cr}, impl::DiscreteBase<V> {rest...},
                   impl::OneDrawBase {} {}

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        void begin_update(const std::size_t idx, const std::size_t cps) override {
            sample_rademacher(cps);
            Parent::begin_update(idx, cps);
        }

        void clear() override {
            Parent::clear();
            u_rade.clear();
        }

        void end_update() override {
            Parent::end_update();
            u_rade.clear();
        }

        /** Merge
         */
        ThisClass& merge(ThisClass& oth) {
            impl::DiscreteBase<V>::do_merge(oth);
            Parent::do_merge(oth);
            return *this;
        }

}; // DISCRETE | ONEDRAW

template <typename T,typename V,typename Hash,typename EqualTo,typename LessV>
class Projection<T,Wheelwright::Features::ONEDRAW |
Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV> :
public impl::ProjectionBase<T,Wheelwright::Features::ONEDRAW |
Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>, public impl::OneDrawBase {
    private:
        using Vector = impl::Vector<V,Wheelwright::Features::ONEDRAW |
            Wheelwright::Features::VANILLA,LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Easy alias for the ProjectionBase parent */
        using Parent = impl::ProjectionBase<T,Wheelwright::Features::ONEDRAW |
            Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>;
        using ThisClass = Projection<T,Wheelwright::Features::ONEDRAW |
            Wheelwright::Features::VANILLA,V,Hash,EqualTo,LessV>;

        void add_zeroes(vec_p &vec) override {
            vec->insert(0, this->u_idx, this->u_cps, u_rade, this->u_avg);
        }

        vec_p built_from_scratch(const V v) const override {
            return std::make_shared<Vector>(v, this->u_idx, this->u_cps, u_rade,
                    v - this->u_avg);
        }

        vec_p extended(const Vector &vec, V v) const override {
            return std::make_shared<Vector>(vec, v, this->u_idx, this->u_cps,
                    u_rade, v - this->u_avg);
        }

    public:
        /** Constructor.
         *
         * @tparam S  variadic template type for the types of the additional
         *            arguments to the constructor.
         *
         * @param cr    the number of references for the catch-all vector. I.e.,
         *              the size of the family of functions, if finite. An
         *              infinite function can be partially simulated by passing
         *              std::numeric_limits<T>::max().
         * @param rest  additional arguments forwarded to the constructore of a
         *              VectorFactory<Vector> object.
         */
        explicit Projection(const std::size_t cr) : Parent {cr},
                 impl::OneDrawBase {} {}

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        void begin_update(const std::size_t idx, const std::size_t cps) override {
            sample_rademacher(cps);
            Parent::begin_update(idx, cps);
        }

        void clear() override {
            Parent::clear();
            u_rade.clear();
        }

        void end_update() override {
            Parent::end_update();
            u_rade.clear();
        }

        /** Merge
         */
        ThisClass& merge(ThisClass& oth) {
            Parent::do_merge(oth);
            return *this;
        }
}; // ONEDRAW | VANILLA

template <typename T,typename V,typename Hash,typename EqualTo, typename LessV>
class Projection<T,Wheelwright::Features::DISCRETE |
Wheelwright::Features::ONEDRAW | Wheelwright::Features::VANILLA,V,Hash,EqualTo,
    LessV> : public impl::ProjectionBase<T,Wheelwright::Features::DISCRETE |
    Wheelwright::Features::ONEDRAW | Wheelwright::Features::VANILLA,V,Hash,
    EqualTo,LessV>, public impl::DiscreteBase<V>, public impl::OneDrawBase {
    private:
        using Vector = impl::Vector<V,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW | Wheelwright::Features::VANILLA,
            LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Easy alias for the ProjectionBase parent */
        using Parent = impl::ProjectionBase<T,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW | Wheelwright::Features::VANILLA,V,
            Hash,EqualTo,LessV>;
        using ThisClass= Projection<T,Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW | Wheelwright::Features::VANILLA,V,
            Hash,EqualTo,LessV>;

        void add_zeroes(vec_p &vec) override {
            vec->insert(0, this->u_idx, this->u_cps, u_rade, this->u_avg);
        }

        vec_p built_from_scratch(const V v) const override {
            return std::make_shared<Vector>(this->codomain.cbegin(),
                    this->codomain.cend(), v, this->u_idx, this->u_cps, u_rade,
                    v - this->u_avg);
        }

        vec_p extended(const Vector &vec, V v) const override {
            return std::make_shared<Vector>(vec, v, this->u_idx, this->u_cps,
                    u_rade, v - this->u_avg);
        }

    public:
        /** Constructor.
         *
         * @tparam S  variadic template type for the types of the additional
         *            arguments to the constructor.
         *
         * @param cr    the number of references for the catch-all vector. I.e.,
         *              the size of the family of functions, if finite. An
         *              infinite function can be partially simulated by passing
         *              std::numeric_limits<T>::max().
         * @param rest  additional arguments forwarded to the constructore of a
         *              VectorFactory<Vector> object.
         */
        template<typename... S> explicit Projection(const std::size_t cr,
                S... rest) : Parent {cr}, impl::DiscreteBase<V> {rest...},
                   impl::OneDrawBase {} {}

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        void begin_update(const std::size_t idx, const std::size_t cps) override {
            sample_rademacher(cps);
            Parent::begin_update(idx, cps);
        }

        void clear() override {
            Parent::clear();
            u_rade.clear();
        }

        void end_update() override {
            Parent::end_update();
            u_rade.clear();
        }

        /** Merge
         */
        ThisClass& merge(ThisClass& oth) {
            impl::DiscreteBase<V>::do_merge(oth);
            Parent::do_merge(oth);
            return *this;
        }
}; // DISCRETE | ONEDRAW | VANILLA

} // namespace Wheelwright

#endif // WHEELWRIGHT_PROJECTION_H
