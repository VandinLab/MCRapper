/** @file Vector.h
 * Projection vector
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

#ifndef WHEELWRIGHT_VECTOR_H
#define WHEELWRIGHT_VECTOR_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "matteoutils/KeyValueGetters.h"
#include "matteoutils/NumericTools.h"

#include "Features.h"

namespace Wheelwright::impl {

enum Rademacher: signed char;

/** Empty template definition for the Vector's in a Projection.
 *
 * The elements of a Projection are vectors with a number of components equal to
 * the current size of the sample for which the Projection is considered. In
 * practice though, the individual components of the vectors are never needed:
 * the sum of the components is needed to compute the sample average, and some
 * other aggregate quantities are used to compute upper bounds to the Rademacher
 * averages or to the empirical Radeacher averages. The aggregate quantities to
 * keep track of are specified using @a F. See the docs for Wheelwright::Features
 * about the possible aggregates.
 *
 * The components of the Vectors are zero-indexed.
 *
 * @tparam T        the type for the components in the Vector.
 * @tparam F        the Wheelwright::Features supported by the Vector.
 * @tparam Less     the less-than functor for @a T. Default: std::less<T>.
 * @tparam J        used for std::enable_if_t tests in specializations.
 */
template<typename T, int F, typename Less = std::less<T>, typename J = T>
class Vector  {
    Vector() = 0;
};

/** Specialization for zero features, acting as abstract base class for all the
 * others.
 *
 *  The sum of the entries of the Vector is maintained in this class.
 *
 *  An hash value is maintained throughout the lifetime of the object to
 *  uniquely (up to collisions) identify a Vector,
 *
 * In order to support progressive sampling algorithms, the vector may grow over
 * time (see insert()), in the sense that additional components can be added.
 *
 * @tparam T        the type for the components in the Vector.
 * @tparam Less     the less-than functor for @a T. Default: std::less<T>.
 */
template<typename T, typename Less> class Vector<T,Wheelwright::Features::ZERO,Less,T>  {
    public:
        using value_type = T;

        /** Returns a hash representing the vector. */
        std::size_t hash() const { return hash_; }

        /** Less-than comparison function.
         *
         * @param e  the Vector to compare to.
         *
         * @returns  true if the hash of @a *this is less than the hash of
         *           @a e.
         */
        bool operator<(const
                Vector<T,Wheelwright::Features::ZERO,Less,T> &e) const {
            return hash_ < e.hash_;
        }

        /** Equality comparison function.
         *
         * It uses the operator< comparison to determine the result.
         *
         * @param e  the Vector to compare to.
         */
        bool operator==(const
                Vector<T,Wheelwright::Features::ZERO,Less,T> &e)
            const {
            return ! (*this < e || e < *this);
        }

        /** Returns the sum of the components in the vector. */
        T sum() const { return sum_; }

    protected:
        /** Constructs a new Vector with the components from @a i to @a i + @a c
         * - 1 set to @a v.
         *
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         */
        Vector(const T v, const std::size_t i, const std::size_t c) :
                sum_ {static_cast<T>(v * c)} {
#ifdef WHEELWRIGHT_INVARIANTS
            if (c == 0)
                throw std::out_of_range("Wheelwright::Vector::Vector: c must "
                        "be positive.");
#endif
            boost::hash_combine(hash_, i);
            boost::hash_combine(hash_, c);
            boost::hash_combine(hash_, v);
        }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        Vector(const Vector<T,Wheelwright::Features::ZERO,Less,T>
                &a,
                const Vector<T,Wheelwright::Features::ZERO,Less,T>
                &b) :
                sum_{a.sum_ + b.sum_} {
            // We want merging of Projections to be symmetric, so we
            // combine the hashes of the Vectors always in the same way.
            if (a.hash_ <= b.hash_) {
                hash_ = a.hash_;
                boost::hash_combine(hash_, b.hash_);
            } else {
                hash_ = b.hash_;
                boost::hash_combine(hash_, a.hash_);
            }
        }

        /** Protected and non-virtual destructor because this class only acts
         * as a base class.
         */
        ~Vector() = default;

        /** Inserts @a c components with value @a v starting from the component
         * with index @a i (included).
         *
         * There is no safeguard preventing the user from adding the same
         * component(s) multiple times, even with different values, or any
         * constraint on the order of the updates. Both aspects have an impact
         * on the hash, and therefore also on the comparison offered by
         * operator<().
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         *
         * @returns the new hash of the Vector.
         */
        virtual std::size_t insert(const T v, const std::size_t i,
                const std::size_t c) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (c == 0)
                throw std::out_of_range("Wheelwright::Vector::insert: c must "
                        "be positive.");
#endif
            boost::hash_combine(hash_, i);
            boost::hash_combine(hash_, c);
            boost::hash_combine(hash_, v);
            sum_ += v * c;
            return hash_;
        }

        // insert() with 4 arguments so it can be overridden, and if the
        // overriding functions set a default value for the last parameter, the
        // compiler does not think that the insert() with 3 arguments is being
        // called.
        virtual std::size_t insert(const T v, const std::size_t i,
                const std::size_t c, const double) {
            return insert(v, i, c);
        }

        /** The hash of the vector. */
        std::size_t hash_ {0};
        /** The sum of the components. */
        T sum_;
};

/** A specialization for Vector offering the minimal functionality to be used
 * as element in a Projection while imposing no restrictions on the domain of
 * the components.
 *
 * The components in the Vector can be from any domain, including continuous
 * ones, i.e. do not need to be from a specific final set (for such
 * restrictions, see the specializations having the template parameter @a F with
 * the bits for Wheelwright::Features::DISCRETE set.)
 *
 * @tparam T        the type for the components in the Vector.
 * @tparam Less     the less-than functor for @a T. Default: std::less<T>.
 */
template<typename T, typename Less> class Vector<T,
    Wheelwright::Features::VANILLA,Less,T> :
    public Vector<T,Wheelwright::Features::ZERO,Less,T> {

    public:
        /** Constructs a new Vector with the components from @a i to @a i + @a c
         * - 1 set to @a v.
         *
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::limits<double>::max(), @a v is used.
         **/
        Vector(const T v, const std::size_t i, const std::size_t c,
                const double lv = std::numeric_limits<double>::max()) :
                parent(v, i, c), elltwosquared_ {(lv ==
                        std::numeric_limits<double>::max()) ?
                    std::pow(v, 2) * c : std::pow(lv, 2) * c} {}

        /** Constructor from an existing Vector<T> with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::numeric_limits<double>::max(), @a v
         *           is used.
         */
        Vector(const Vector<T, Wheelwright::Features::VANILLA> &e,
                const T v, const std::size_t i, const std::size_t c,
                const double lv = std::numeric_limits<double>::max())
            : Vector(e) {
            insert(v, i, c, lv);
        }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        Vector(const Vector<T,Wheelwright::Features::VANILLA,Less,T> &a,
                const Vector<T,Wheelwright::Features::VANILLA,Less,T> &b) :
            parent(a, b), elltwosquared_ {a.elltwosquared_ + b.elltwosquared_} {}

        /** Virtual destructor for base class. */
        virtual ~Vector() = default;

        /** Returns the \f$\ell_2\f$-norm of the vector. */
        double elltwo() const { return std::sqrt(elltwosquared_); }

        /** Returns the squared \f$\ell_2\f$-norm of the vector. */
        double elltwosquared() const { return elltwosquared_; }

        /** Inserts @a c components with value @a v starting from the component
         * with index @a i (included).
         *
         * There is no safeguard preventing the user from adding the same
         * component(s) multiple times, even with different values, or any
         * constraint on the order of the updates. Both aspects have an impact
         * on the hash, and therefore also on the comparison offered by
         * operator<().
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::numeric_limits<double>::max(), @a v
         *           is used.
         *
         * @returns the new hash of the Vector.
         */
        virtual std::size_t insert(const T v, const std::size_t i,
                const std::size_t c,
                const double lv = std::numeric_limits<double>::max()) override {
            if (lv == std::numeric_limits<double>::max())
                elltwosquared_ += std::pow(v, 2) * c;
            else
                elltwosquared_ += std::pow(lv, 2) * c;
            return parent::insert(v, i, c);
        }

    private:
        /** Base class alias for easier code reading/writing. */
        using parent = Vector<T,Wheelwright::Features::ZERO,Less,T>;
        /** The squared ell-2 norm of the vector. */
        double elltwosquared_;
}; // class Vector

/** A second abstract base class for Vector specialization with the template
 * parameter @a F with the bits of Wheelwright::Features::DISCRETE set.
 *
 * The main purpose of this class is to keep track, for each value in the
 * discrete finite domain of the components, how many components have that
 * value.
 *
 * @tparam T        the type for the components in the Vector.
 * @tparam Less     the less-than functor for @a T. Default: std::less<T>.
 */
template<typename T, typename Less = std::less<T>> class DiscreteVectorBase {

    public:
        /** Return an iterator to the first count. */
        auto begin() const { return counts_begin(); }

        /** Returns the number of components with value @a v.
         *
         * The function does not check whether @a v is in the domain of
         * the components.
         *
         * @param v  the value whose count to return.
         */
        std::size_t count(const T &v) const { return vc.at(v); }

        /** Return an iterator to the first count. */
        auto counts_begin() const {
            return boost::make_transform_iterator(vc.cbegin(),
                    MatteoUtils::ValueGetter<decltype(vc.cbegin())>());
        }

        /** Returns an iterator past the last count. */
        auto counts_end() const {
            return boost::make_transform_iterator(vc.cend(),
                    MatteoUtils::ValueGetter<decltype(vc.cend())>());
        }

        /** Returns an iterator to the first element of the domain.
         *
         * The order of the elements of the domain is determined by @a Less.
         */
        auto domain_begin() const {
            return boost::make_transform_iterator(vc.cbegin(),
                    MatteoUtils::KeyGetter<decltype(vc.cbegin())>());
        }

        /** Returns an iterator past the last element of the domain.
         *
         * The order of the elements of the domain is determined by @a Less.
         */
        auto domain_end() const {
            return boost::make_transform_iterator(vc.cend(),
                    MatteoUtils::KeyGetter<decltype(vc.cend())>());
        }

        /** Returns an iterator past the last count. */
        auto end() const { return counts_end(); }

    protected:
        /** Constructs a new DiscreteVectorBase whose components will belong to
         * the domain specified in the sequence [@ab, @e), and with @a c initial
         * components set to @a v.
         *
         * For the purpose of this class, the indices of the initial components
         * do not matter.
         *
         * @param v  the value for the initial components.
         * @param c  the number of initial components.
         */
        template<typename Iterator> DiscreteVectorBase(Iterator b,
                Iterator e, const T &v, const std::size_t c) {
            static_assert(std::is_same<
                    typename std::iterator_traits<Iterator>::value_type,
                    T>::value, "Wheelwright::impl::DiscreteVectorBase::"
                    "DiscreteVectorBase: the value_type of the type Iterator "
                    "for b and e is not the same as the type T for v.");
            bool found {false};
            for (; b != e; ++b) {
                if (*b != 0) {
                    std::size_t cnt {0};
                    if (! found && MatteoUtils::logically_equal(v, *b, 2)) {
                        found = true;
                        cnt = c;
                    }
#ifdef WHEELWRIGHT_INVARIANTS
                    if (! vc.emplace(*b, cnt).second) {
                        std::stringstream ss;
                        ss << "Wheelwright::impl::DiscreteVectorBase::"
                            "DiscreteVectorBase: duplicate element '" << *b
                            << "' in the passed sequence for the components "
                            "domain.";
                        throw std::runtime_error(ss.str());
                    }
#else
                    vc.emplace(*b, cnt).second;
#endif
                }
            }
#ifdef WHEELWRIGHT_INVARIANTS
            if (vc.empty())
                throw std::runtime_error(
                        "Wheelwright::impl::DiscreteVectorBase::"
                        "DiscreteVectorBase: the passed sequence for the "
                        "components domain does not contain non-zero "
                        "elements.");
            if (! found && v != 0) {
                std::stringstream ss;
                ss << "Wheelwright::impl::DiscreteVectorBase::"
                    "DiscreteVectorBase: v=" << v << " is not in the "
                    "components domain {0, ";
                for (const auto &p : vc)
                    ss << p.first << ", ";
                ss.seekp(-2, ss.cur);
                ss << "}.";
                throw std::runtime_error(ss.str());
            }
#endif
        }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        DiscreteVectorBase(const DiscreteVectorBase<T,Less> &a,
                const DiscreteVectorBase<T,Less> &b) : vc {a.vc.begin(),
            a.vc.end()} {
#ifdef WHEELWRIGHT_INVARIANTS
            if (vc.size() != b.vc.size())
                throw std::runtime_error("Wheelwright::DiscrVector::"
                        "DiscreteVectorBase: the two vectors have different "
                        "maps.");
#endif
            for (auto&& valco : vc) {
                auto it {b.vc.find(valco.first)};
#ifdef WHEELWRIGHT_INVARIANTS
                if (it == b.vc.end())
                    throw std::runtime_error("Wheelwright::DiscrVector::"
                        "DiscreteVectorBase: the two vectors have different "
                        "maps.");
#endif
                 valco.second += it->second;
            }
        }

        /** Inserts @a c components with value @a v.
         *
         * For the purpose of this class, the indices of the added components do
         * not matter.
         *
         * @param v  the value for the additional components.
         * @param c  the number of additional components.
         */
        void insert(const T &v, const std::size_t c) {
            if (auto it {vc.find(v)}; it == vc.end()) {
                for (it = vc.begin(); it != vc.end() &&
                        ! MatteoUtils::logically_equal(v, it->first, 2); ++it)
                    ; // empty-body for loop
                if (it == vc.end() && v != 0) {
                    std::stringstream ss;
                    ss << "Wheelwright::impl::DiscreteVectorBase::insert: v="
                        << v << " is not in the components domain {0, ";
                    for (const auto &p : vc)
                        ss << p.first << ", ";
                    ss.seekp(-2, ss.cur);
                    ss << "}.";
                    throw std::runtime_error(ss.str());
                }
            } else
                it->second += c;
        }

        /** Protected and non-virtual destructor because this class only acts as
         * base class.
         */
        ~DiscreteVectorBase() = default;

    private:
        // XXX Wheelwright::impl::DiscreteVectorBase: era_bound_AGOR_FD relies
        // on the fact that the vc container is "sorted" (i.e., not a
        // unordered_map). (20170911: I don't remember why there is this
        // assumption in era_bound_AGOR_FD.)
        /** The map containing, for each value of the domain, the number of
         * components with that value.
         */
        std::map<T, std::size_t,Less> vc;
};

/** A specialization for Vector offering the minimal functionality to be used
 * as element in a Projection while imposing restrictions on the domain of
 * the components.
 *
 * The components in the Vector must be from a finite discrete domain specified
 * at construction time (for a class without such restrictions, see the
 * specializations having the template parameter @a F with the bits for
 * Wheelwright::Features::VANILLA set.)
 *
 * @tparam T        the type for the components in the Vector.
 * @tparam Less     the less-than functor for @a T. Default: std::less<T>.
 */
template<typename T, int F, typename Less> class Vector<T,F,Less,
    std::enable_if_t<F == Wheelwright::Features::DISCRETE, T>> :
    public Vector<T, F ^ Wheelwright::Features::DISCRETE,Less,T>,
    public DiscreteVectorBase<T,Less> {

    private:
        using DiscreteVectorBase<T,Less>::insert;
        /** Base class alias for easier code reading/writing. */
        using parent = Vector<T, F ^ Wheelwright::Features::DISCRETE,
                   Less,T>;

    public:
        /** Constructor explicitly specifying the non-zero elements of the
         * domain of the components and the initial components.
         *
         * Constructs a vector whose components belong to the finite discrete
         * domain specified in the sequence [@a b, @a e), with initial
         * components given in the remaining arguments.
         *
         * @tparam Iterator  the Iterator type for the arguments @a b and @a e.
         *                    The value_type trait must be T.
         *
         * @param [b,e)  the sequence specifying the non-zero elements of the
         *               components domain. An eventual zero is ignored.
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         */
        template<typename Iterator> Vector(Iterator b, Iterator e, const T &v,
                const std::size_t i, const std::size_t c) : parent(v, i, c),
        DiscreteVectorBase<T,Less>(b, e, v, c) { }

        /** Constructor from an existing Vector with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         */
        Vector(const Vector<T,F,Less,T> &e, const T &v,
                const std::size_t i, const std::size_t c) : Vector(e) {
            insert(v, i, c);
        }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        Vector(const Vector<T,F,Less,T> &a, const Vector<T,F,Less,T> &b) :
                Vector<T,
                    F ^ Wheelwright::Features::DISCRETE,Less,
                    T>(a,b), DiscreteVectorBase<T,Less>(a,b) { }

        /** Virtual destructor for a base class. */
        virtual ~Vector() = default;

        /** Inserts components into the vector.
         *
         * There is no safeguard preventing the user from adding the same
         * component(s) multiple times, even with different values, or any
         * constraint on the order of the updates. Both aspects have an impact
         * on the hash, and therefore also on the comparison offered by
         * operator<().
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         *
         * @returns  the new hash of the vector.
         */
        virtual std::size_t insert(const T v, const std::size_t i,
                const std::size_t c) override {
            insert(v, c);
            return parent::insert(v, i, c);
        }
};

/** A specialization for Vector offering the functionalities of both VANILLA
 * Vectors and DISCRETE vectors.
 *
 * The components in the Vector must be from a finite discrete domain specified
 * at construction time, but the @a insert() member function has the same
 * arguments as the one for VANILLA Vectors, allowing to derive different kinds
 * of bounds to the empirical Rademacher averages.
 *
 * @tparam T        the type for the components in the Vector.
 * @tparam Less     the less-than functor for @a T. Default: std::less<T>.
 */
template<typename T, int F, typename Less> class Vector<T,F,Less,
    std::enable_if_t<F == (Wheelwright::Features::DISCRETE |
            Wheelwright::Features::VANILLA), T>> : public
    Vector<T, F ^ Wheelwright::Features::DISCRETE,Less,T>,
    public DiscreteVectorBase<T,Less> {

    private:
        using DiscreteVectorBase<T,Less>::insert;
        /** Base class alias for easier code reading/writing. */
        using parent = Vector<T, F ^ Wheelwright::Features::DISCRETE,
                       Less>;

    public:
        /** Constructor explicitly specifying the non-zero elements of the
         * domain of the components and the initial components.
         *
         * Constructs a vector whose components belong to the finite discrete
         * domain specified in the sequence [@a b, @a e), with initial
         * components given in the remaining arguments.
         *
         * @tparam Iterator  the Iterator type for the arguments @a b and @a e.
         *                    The value_type trait must be T.
         *
         * @param [b,e)  the sequence specifying the non-zero elements of the
         *               components domain. An eventual zero is ignored.
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::numeric_limits<double>::max(), @a v
         *           is used.
         */
        template<typename Iterator> Vector(Iterator b, Iterator e, const T &v,
                const std::size_t i, const std::size_t c,
                const double lv = std::numeric_limits<double>::max()) :
            parent(v, i, c, lv), DiscreteVectorBase<T,Less>(b, e, v, c) { }

        /** Constructor from an existing Vector with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         */
        Vector(const Vector<T,F,Less,T> &e, const T &v, const std::size_t i,
                const std::size_t c,
                const double lv = std::numeric_limits<double>::max())
            : Vector<T,F,Less,T>(e) {
            insert(v, i, c, lv);
        }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        Vector(const Vector<T,F,Less,T> &a, const Vector<T,F,Less,T> &b) :
                Vector<T,
                    F ^ Wheelwright::Features::DISCRETE,Less,
                    T>(a,b), DiscreteVectorBase<T,Less>(a,b) { }

        /** Inserts @a c components with value @a v starting from the component
         * with index @a i (included).
         *
         * There is no safeguard preventing the user from adding the same
         * component(s) multiple times, even with different values, or any
         * constraint on the order of the updates. Both aspects have an impact
         * on the hash, and therefore also on the comparison offered by
         * operator<().
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::numeric_limits<double>::max(), @a v
         *           is used.
         *
         * @returns the new hash of the Vector.
         */
        std::size_t insert(const T v, const std::size_t i,
                const std::size_t c,
                const double lv = std::numeric_limits<double>::max()) override final {
            insert(v, c);
            return parent::insert(v, i, c, lv);
        }
};

/** A second abstract base class for Vector specialization with the template
 * parameter @a F with the bits of Wheelwright::Features::ONEDRAW set.
 *
 * The main purpose of this class is to keep track of the sum of the products of
 * the components with the Rademacher variables.
 *
 * @tparam T        the type for the components in the Vector.
 */
template<typename T> class OneDrawVectorBase {
    public:
        /** Returns the sum of the products of the components with the
         * Rademacher variables.
         */
        double onedrawsum() const { return onedrawsum_; }

    protected:

        OneDrawVectorBase(const T v, const std::size_t c,
                const std::vector<Rademacher> &r) {
            insert(v, c, r);
        }

        OneDrawVectorBase(const double s = 0) : onedrawsum_ {s} {}

        /** Protected and non-virtual destructor because this class only acts as
         * a base class.
         */
        ~OneDrawVectorBase() = default;

        /** Inserts @a c components with value @a v, multiplied by the
         * Rademacher variables in @a r.
         *
         * For the purpose of this class, the indices of the added components do
         * not matter.
         *
         * @param v  the value for the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        void insert(const T v, const std::size_t c,
                const std::vector<Rademacher> &r) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (c != r.size())
                throw std::invalid_argument(
                        "Wheelwright::impl::OneDrawVectorBase::insert: the "
                        "length of r must be equal to c.");
            if (c == 0)
                throw std::invalid_argument(
                        "Wheelwright::impl::OneDrawVectorBase::insert: the "
                        "number of copies cannot be zero.");
#endif
            for (std::size_t i {0}; i < c; ++i)
                onedrawsum_ += static_cast<double>(v) * r.at(i);
        }

    private:
        /** The sum of the components weighted by the Rademacher averages. */
        double onedrawsum_ {0};
};

/** TODO: Wheelwright::Vector: document ONEDRAW Vectors. */
template<typename T, int F, typename Less> class Vector<T,F,Less,
      std::enable_if_t<F == Wheelwright::Features::ONEDRAW,T>> :
      public Vector<T, F ^ Wheelwright::Features::ONEDRAW>,
      public OneDrawVectorBase<T> {
    private:
        /** Base class alias for easier code reading/writing. */
        using parent = Vector<T,
                    F ^ Wheelwright::Features::ONEDRAW,Less,
                    T>;
        // Needed to avoid warning of hiding an overloaded virtual function
        using parent::insert;

    public:
        /** Constructs new vector with initial components.
         *
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         * @param r  a std::vector of Rademacher variables.
         */
        Vector(const T &v, const std::size_t i, const std::size_t c,
                const std::vector<Rademacher> &r) :
            parent {v, i, c}, OneDrawVectorBase<T> {v, c, r} {
#ifdef WHEELWRIGHT_INVARIANTS
                if (c == 0)
                    throw std::out_of_range(
                            "Wheelwright::impl::Vector::Vector: c must be "
                            "positive.");
#endif
        }

        /** Constructor from an existing Vector with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        Vector(const Vector<T,F,Less,T> &e, const T
                &v, const std::size_t i,
                const std::size_t c,
                const std::vector<Rademacher> &r) : Vector(e) {
#ifdef WHEELWRIGHT_INVARIANTS
                if (c == 0)
                    throw std::out_of_range(
                            "Wheelwright::impl::Vector::Vector: c must be "
                            "positive.");
#endif
                insert(v, i, c, r);
            }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        Vector(const Vector<T,F,Less,T> &a, const Vector<T,F,Less,T> &b) :
                parent(a, b),
                OneDrawVectorBase<T>(a.onedrawsum() + b.onedrawsum()) {}

        /** Virtual destructor for a base class. */
        virtual ~Vector() = default;

        /** Inserts @a c components with value @a v starting from the component
         * with index @a i (included).
         *
         * For the purpose of this class, the indices of the added components do
         * not matter.
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        virtual std::size_t insert(const T v, const std::size_t i,
                const std::size_t c,
                const std::vector<Rademacher> &r) {
            OneDrawVectorBase<T>::insert(v, c, r);
            return parent::insert(v, i, c);
        }
};

/** A specialization for Vector offering the functionalities of DISCRETE Vectors
 * and ONEDRAW Vectors.
 *
 * TODO: Wheelwright::Vector: Complete the documentation for DISCRETE|ONEDRAW
 */
template<typename T, int F, typename Less> class Vector<T,F,Less,
    std::enable_if_t<F == (Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW), T>> :
    public Vector<T, F ^ Wheelwright::Features::DISCRETE,Less,T>,
    public DiscreteVectorBase<T,Less> {

    private:
        /** Base class alias for easier code reading/writing. */
        using parent = Vector<T, F ^ Wheelwright::Features::DISCRETE,Less>;

    public:
        /** Constructor explicitly specifying the non-zero elements of the
         * domain of the components and the initial components.
         *
         * Constructs a vector whose components belong to the finite discrete
         * domain specified in the sequence [@a b, @a e), with initial
         * components given in the remaining arguments.
         *
         * @tparam Iterator  the Iterator type for the arguments @a b and @a e.
         *                    The value_type trait must be T.
         *
         * @param [b,e)  the sequence specifying the non-zero elements of the
         *               components domain. An eventual zero is ignored.
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         * @param r  a std::vector of Rademacher variables.
         */
        template<typename Iterator> Vector(Iterator b, Iterator e,
                const T v, const std::size_t i, const std::size_t c,
                const std::vector<Rademacher> &r) :
            parent(v, i, c, r), DiscreteVectorBase<T,Less>(b, e, v, c) {}

        /** Constructor from an existing Vector with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        Vector(const Vector<T,F,Less> &e, const T &v, const std::size_t i,
                const std::size_t c,
                const std::vector<Rademacher> &r) : Vector(e) {
            insert(v, i, c, r);
        }

        /** Constructor from an existing Vector with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         */
        Vector(const Vector<T,F,Less,T> &a, const Vector<T,F,Less,T> &b) :
                parent(a, b), DiscreteVectorBase<T,Less>(a, b) { }

        /** Virtual destructor for a base class. */
        virtual ~Vector() = default;

        /** Inserts @a c components with value @a v starting from the component
         * with index @a i (included).
         *
         * For the purpose of this class, the indices of the added components do
         * not matter.
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        std::size_t insert(const T v, const std::size_t i,
                const std::size_t c,
                const std::vector<Rademacher> &r) override final {
            DiscreteVectorBase<T,Less>::insert(v, c);
            return parent::insert(v, i, c, r);
        }
};

/** A specialization for Vector offering the functionalities of VANILLA
 * Vectors and ONEDRAW Vectors.
 *
 * TODO: Wheelwright::Vector: Complete the documentation for VANILLA|ONEDRAW
 */
template<typename T, int F, typename Less> class Vector<T,F,Less,
      std::enable_if_t<F == (Wheelwright::Features::ONEDRAW |
              Wheelwright::Features::VANILLA),T>> :
      public Vector<T, F ^ Wheelwright::Features::ONEDRAW>,
      public OneDrawVectorBase<T> {

    private:
        /** Base class alias for easier code reading/writing. */
        using parent = Vector<T, F ^ Wheelwright::Features::ONEDRAW,Less,T>;
        // Needed to avoid warning of hiding an overloaded virtual function
        using parent::insert;

    public:
        /** Constructs a new Vector with initial components.
         *
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         * @param r  a std::vector of Rademacher variables.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::limits<double>::max(), @a v is used.
         */
        Vector(const T &v, const std::size_t i, const std::size_t c,
                const std::vector<Rademacher> &r,
                const double lv = std::numeric_limits<double>::max()) :
            parent {v, i, c, lv}, OneDrawVectorBase<T> {v, c, r} {
#ifdef WHEELWRIGHT_INVARIANTS
                if (c == 0)
                    throw std::out_of_range(
                            "Wheelwright::impl::Vector::Vector: c must be "
                            "positive.");
#endif
        }
        /** Constructor from an existing Vector with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::limits<double>::max(), @a v is used.
         */
        Vector(const Vector<T,F,Less,T> &e, const T v, const
                std::size_t i, const std::size_t c,
                const std::vector<Rademacher> &r,
                const double lv = std::numeric_limits<double>::max())
            : Vector<T,F,Less,T>(e) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (c == 0)
                throw std::out_of_range(
                        "Wheelwright::impl::Vector::Vector: c must be "
                        "positive.");
#endif
            insert(v, i, c, r, lv);
        }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        Vector(const Vector<T,F,Less,T> &a, const Vector<T,F,Less,T> &b) :
                parent(a, b),
                OneDrawVectorBase<T> {a.onedrawsum() + b.onedrawsum()} {}

        /** Virtual destructor for a base class. */
        virtual ~Vector() = default;

        /** Inserts @a c components with value @a v starting from the component
         * with index @a i (included).
         *
         * For the purpose of this class, the indices of the added components do
         * not matter.
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        virtual std::size_t insert(const T v,
                const std::size_t i, const std::size_t c,
                const std::vector<Rademacher> &r,
                const double lv = std::numeric_limits<double>::max()) {
            OneDrawVectorBase<T>::insert(v, c, r);
            return parent::insert(v, i, c, lv);
        }
};

/** A specialization for Vector offering the functionalities of VANILLA
 * Vectors, ONEDRAW Vectors, and DISCRETE Vectors.
 *
 * TODO: Wheelwright::Vector: Complete the documentation for VANILLA|ONEDRAW|DISCRETE
 */
template<typename T, int F, typename Less> class Vector<T,F,Less,
    std::enable_if_t<F == (Wheelwright::Features::DISCRETE |
            Wheelwright::Features::ONEDRAW |
            Wheelwright::Features::VANILLA), T>> :
    public Vector<T, F ^ Wheelwright::Features::DISCRETE,Less,T>,
    public DiscreteVectorBase<T,Less> {

    private:
        /** Base class alias for easier code reading/writing. */
        using parent = Vector<T, F ^ Wheelwright::Features::DISCRETE,Less,T>;

    public:
        /** Constructor explicitly specifying the non-zero elements of the
         * domain of the components and the initial components.
         *
         * Constructs a vector whose components belong to the finite discrete
         * domain specified in the sequence [@a b, @a e), with initial
         * components given in the remaining arguments.
         *
         * @tparam Iterator  the Iterator type for the arguments @a b and @a e.
         *                    The value_type trait must be T.
         *
         * @param [b,e)  the sequence specifying the non-zero elements of the
         *               components domain. An eventual zero is ignored.
         * @param v  the value for the initial components.
         * @param i  the starting index of the initial components.
         * @param c  the number of initial components.
         * @param r  a std::vector of Rademacher variables.
         * @param lv the "centralized" (non-squared) value to use for the
         *           elltwosquared. If std::limits<double>::max(), @a v is used.
         */
        template<typename Iterator> Vector(Iterator b, Iterator e, const T v,
                const std::size_t i, const std::size_t c,
                const std::vector<Rademacher> &r,
                const double lv = std::numeric_limits<double>::max()) :
            parent(v, i, c, r, lv), DiscreteVectorBase<T,Less>(b, e, v, c) {}

        /** Constructor from an existing Vector with some additional
         * components.
         *
         * Constructs a vector which has the same components as @a e, plus some
         * additional components specified in the arguments.
         *
         * @param e  the "base" vector.
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        Vector(const Vector<T,F,Less,T> &e, const T v, const std::size_t i,
                const std::size_t c,
                const std::vector<Rademacher> &r,
                const double lv = std::numeric_limits<double>::max()) : Vector(e) {
            insert(v, i, c, r, lv);
        }

        /** Constructs a new Vector as the concatenation of the Vector @a a and
         * the Vector @a b.
         *
         * The order of the arguments is irrelevant: a call with (@a a, @a b)
         * constructs a Vector equivalent to a call with (@a b, @a a).
         *
         * @param a,b  the Vectors to merge.
         */
        Vector(const Vector<T,F,Less,T> &a, const Vector<T,F,Less,T> &b) :
                parent(a,b), DiscreteVectorBase<T,Less>(a,b) { }

        /** Inserts @a c components with value @a v starting from the component
         * with index @a i (included).
         *
         * For the purpose of this class, the indices of the added components do
         * not matter.
         *
         * @param v  the value for the additional components.
         * @param i  the starting index of the additional components.
         * @param c  the number of additional components.
         * @param r  a std::vector of Rademacher variables.
         */
        std::size_t insert(const T v, const std::size_t i,
                const std::size_t c,
                const std::vector<Rademacher> &r,
                const double lv) override final {
            DiscreteVectorBase<T,Less>::insert(v, c);
            return parent::insert(v, i, c, r, lv);
        }
};

} // namespace Wheelwright::impl


namespace std {

    /** Specialization of std::hash for Vector.
     *
     * @tparam T        the type for the components in the Vector.
     * @tparam F        the Wheelwright::Features supported by the Vector.
     * @tparam Less     the less-than functor for @a T.
     */
    template<typename T,int F, typename Less>
        struct hash<Wheelwright::impl::Vector<T,F,Less,T>> {
        /** The type of the argument of operator(). */
        using argument_type = Wheelwright::impl::Vector<T,F,Less,T>;
        /** The return type of operator(). */
        using result_type = std::size_t;
        /** Returns the hash of the Vector @a s.
         *
         * @param s  the Vector whose object to return the hash of.
         */
        result_type operator()(argument_type const& s) const {
            return s.hash();
        }
    };
} // namespace std

#endif // WHEELWRIGHT_VECTOR_H
