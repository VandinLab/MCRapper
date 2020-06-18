/** @file Projection_impl.h
 * Projection set.
 *
 * @author Matteo Riondato
 * @date 2018 03 14
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

#ifndef WHEELWRIGHT_PROJECTION_IMPL_H
#define WHEELWRIGHT_PROJECTION_IMPL_H

#include <cstddef>
#include <functional>
#include <memory>

#include "ProjectionIterator.h"
#include "Vector.h"

namespace Wheelwright::impl {

/** A Rademacher random variable. */
enum Rademacher: signed char {
    /** Positive value i.e., 1. */
    POS = 1,
    /** Negative value, i.e., -1. */
    NEG = -1
};


template <typename T, int F, typename V = double, typename Hash = std::hash<T>,
         typename EqualTo = std::equal_to<T>, typename LessV = std::less<V>>
             class ProjectionBase {
    private:
        using Vector = impl::Vector<V,F,LessV>;
        /** Type for pointers to Vector. */
        using vec_p = std::shared_ptr<Vector>;
        /** Type of the const iterator. */
        using const_iterator = impl::ProjectionIterator<typename
            std::unordered_map<vec_p,std::size_t, MatteoUtils::PtrHash<Vector>,
            MatteoUtils::PtrEqualTo<Vector>>::const_iterator>;

        std::size_t dims_ {0};
        std::size_t famsz_;
        std::size_t u_pars {0};
        std::size_t u_vecs {0};
        vec_p cav {nullptr};
        std::map<V, std::vector<T>, LessV> newvals;
        std::unordered_map<vec_p, decltype(newvals),
            MatteoUtils::PtrHash<Vector>,
            MatteoUtils::PtrEqualTo<Vector>> old2vals;
        // TODO: Wheelwright::Projection: instead of using par2vec and vecrefs,
        // consider using Boost.Bimap.
        std::unordered_map<T,vec_p,Hash,EqualTo> par2vec;
        std::unordered_map<vec_p, std::size_t, MatteoUtils::PtrHash<Vector>,
            MatteoUtils::PtrEqualTo<Vector>> u_vecrefs;
        std::unordered_map<vec_p, std::size_t, MatteoUtils::PtrHash<Vector>,
            MatteoUtils::PtrEqualTo<Vector>> vecrefs;

    protected:
        virtual void add_zeroes(vec_p &v) = 0;

        virtual std::shared_ptr<Vector> built_from_scratch(const V v) const = 0;

        virtual std::shared_ptr<Vector> extended(const Vector &vec, V v) const = 0;

        /** Merge a Projection into this one.
         *
         * This function is useful, for example, in OpenMP reductions.
         *
         * If an update is active in @a oth, that update is lost, while an
         * update active in *this is not lost.
         *
         * The @a oth object is invalidated after calling this function and
         * must not be used any further.
         *
         * @param oth  the Projection to merge into this one.
         */
        void do_merge(ProjectionBase<T,F,V,Hash,EqualTo> &oth) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (famsz_ != oth.famsz_)
                throw std::runtime_error("Wheelwright::Projection::do_merge: "
                        "the sizes of the families of functions passed in the "
                        "ctors of the two Projections are not equal.");
#endif
            if (oth.vecrefs.empty()) {
                oth.clear();
                return;
            }
            if (vecrefs.empty()) {
                *this = std::move(oth);
                return;
            }
            std::size_t cav_refs {0};
#ifdef WHEELWRIGHT_INVARIANTS
            std::size_t oth_cav_refs {0};
#endif
            if (cav != nullptr)
                cav_refs = vecrefs.at(cav);
#ifdef WHEELWRIGHT_INVARIANTS
            if (oth.cav != nullptr)
                oth_cav_refs = oth.vecrefs.at(oth.cav);
#endif
            vecrefs.clear();
            for (auto& o2v: par2vec) {
                vec_p vec_n {o2v.second};
                // Merge the vectors that this parameter is mapped to in the two
                // Projections (the one in oth may be the cav of that
                // Projection).
                if (const auto oth_it {oth.par2vec.find(o2v.first)};
                        oth_it != oth.par2vec.end()) {
                    vec_n = std::make_shared<Vector>(*vec_n, *(oth_it->second));
                    // Removing the parameter from oth.par2vec but not
                    // decrementing oth.vecrefs[oth_it->second], thus breaking
                    // the invariants for oth.
                    oth.par2vec.erase(oth_it);
                } else {
#ifdef WHEELWRIGHT_INVARIANTS
                    if (oth_cav_refs == 0)
                        throw std::runtime_error(
                                "Wheelwright::Projection::do_merge: found "
                                "parameter in this Projection that is not "
                                "mapped in the oth Projection, but already "
                                "reached the size of the family of functions "
                                "specified in the call to the ctor of oth.");
                    --oth_cav_refs;
#endif
                    vec_n = std::make_shared<Vector>(*vec_n, *(oth.cav));
                }
                // Insert the vector for this parameter into vecrefs, or just
                // increase its refcount.
                if (auto r {vecrefs.try_emplace(vec_n, 1)}; r.second == true)
                    o2v.second = r.first->first;
                else {
                    // TODO: Wheelwright::Projection::end_update: devise a test
                    // case for this piece of code.
                    r.first->second += 1;
                    o2v.second = r.first->first;
                }
            }
            // Insert into this Projection the parameters that are only present
            // in oth.
#ifdef WHEELWRIGHT_INVARIANTS
            if (cav_refs < oth.par2vec.size())
                throw std::runtime_error(
                        "Wheelwright::Projection::do_merge: the number of "
                        "parameters in the oth Projection that are not mapped "
                        "in this Projection would make the number of "
                        "parameters in this Projection larger than the size of "
                        "the family of " "functions specified in the call to "
                        "the ctor of this.");
#endif
            for (const auto &o2v : oth.par2vec) {
                    // The vector we try to emplace in vecrefs is the merged
                    // vector obtained by merging the vector that this parameter
                    // is mapped to in the oth Projection with the cav of this
                    // Projection.
                if (auto r
                        {vecrefs.try_emplace(
                                std::make_shared<Vector>(*(o2v.second), *(cav)),
                                1)};
                        r.second == true)
                    par2vec.try_emplace(o2v.first, r.first->first);
                else {
                    par2vec.try_emplace(o2v.first, r.first->first);
                    r.first->second += 1;
                }
            }
            cav_refs -= oth.par2vec.size();
#ifdef WHEELWRIGHT_INVARIANTS
            if (cav_refs != oth_cav_refs)
                throw std::runtime_error("Wheelwright::Projection::do_merge: "
                        "the number of references to the catch-all vector do "
                        "not match.");
#endif
            if (cav_refs > 0)
                cav = vecrefs.try_emplace(std::make_shared<Vector>(*cav,
                            *(oth.cav)),
                        cav_refs).first->first;
            else
                cav = nullptr;
            if (oth.dims_ > dims_)
                dims_ = oth.dims_;
            oth.clear();
        }


        std::size_t u_cps {0};
        std::size_t u_idx {0};
        double u_avg {0};

    public:
        /** Constructor.
         *
         * @param cr  the number of references for the catch-all vector.
         *            I.e., the size of the family of functions, if finite. An
         *            infinite function can be partially simulated by passing
         *            std::numeric_limits<T>::max().
         */
        explicit ProjectionBase(const std::size_t cr) : famsz_ {cr} {}

        /** Copy constructor. */
        ProjectionBase(const ProjectionBase<T, F, V, Hash, EqualTo, LessV> &oth) :
                dims_ {oth.dims_}, famsz_ {oth.famsz_}, u_pars {oth.u_pars},
                u_vecs {oth.u_vecs},
                newvals {oth.newvals.begin(), oth.newvals.end()},
                u_cps {oth.u_cps}, u_idx {oth.u_idx}, u_avg {oth.u_avg} {
            std::unordered_map<vec_p, vec_p, MatteoUtils::PtrHash<Vector>,
                MatteoUtils::PtrEqualTo<Vector>> oth2this;
            oth2this.reserve(oth.vecrefs.size());
            vecrefs.reserve(oth.vecrefs.size());
            par2vec.reserve(oth.par2vec.size());
            u_vecrefs.reserve(oth.u_vecrefs.size());
            old2vals.reserve(oth.old2vals.size());
            if (oth.cav != nullptr) {
                cav = vecrefs.try_emplace(std::make_shared<Vector>(*(oth.cav)),
                        oth.vecrefs.at(oth.cav)).first->first;
                oth2this.try_emplace(oth.cav, cav);
            }
            for (const auto &p : oth.vecrefs)
                oth2this.try_emplace(p.first,
                        vecrefs.try_emplace(std::make_shared<Vector>(
                                *(p.first)),
                            p.second).first->first);
            for (const auto &p : oth.par2vec)
                par2vec.try_emplace(p.first, oth2this.at(p.second));
            for (const auto &p : oth.u_vecrefs)
                u_vecrefs.try_emplace(oth2this.at(p.first), p.second);
            for (const auto &p : oth.old2vals)
                old2vals.try_emplace(oth2this.at(p.first), p.second);
        }

        /** Move constructor. */
        ProjectionBase(ProjectionBase<T, F, V, Hash, EqualTo, LessV> &&oth)
            noexcept : dims_ {oth.dims_},
                     famsz_ {oth.famsz_}, u_pars {oth.u_pars},
                     u_vecs {oth.u_vecs}, cav {std::move(oth.cav)},
                     newvals {std::move(oth.newvals)},
                     old2vals {std::move(oth.old2vals)},
                     par2vec {std::move(oth.par2vec)},
                     u_vecrefs {std::move(oth.u_vecrefs)},
                     vecrefs {std::move(oth.vecrefs)}, u_cps {oth.u_cps},
                     u_idx {oth.u_idx}, u_avg {oth.u_avg} { }

        /** Copy assignment operator. */
        ProjectionBase<T, F, V, Hash, EqualTo, LessV>& operator=(
                const ProjectionBase<T, F, V, Hash, EqualTo, LessV> &oth) {
            // Because the class is abstract, we cannot instatiate object from
            // it, so this shortcut of reusing copy-constructor and move
            // assignment won't work.
            // ProjectionBase<T, F, V, Hash, EqualTo, LessV> tmp(oth); // re-use copy-constructor
            //  *this = std::move(tmp); // re-use move-assignment
            dims_  = oth.dims_;
            famsz_  = oth.famsz_;
            u_cps  = oth.u_cps;
            u_idx  = oth.u_idx;
            u_avg  = oth.u_avg;
            u_pars  = oth.u_pars;
            u_vecs  = oth.u_vecs;
            newvals = oth.newvals;
            std::unordered_map<vec_p, vec_p, MatteoUtils::PtrHash<Vector>,
                MatteoUtils::PtrEqualTo<Vector>> oth2this;
            oth2this.reserve(oth.vecrefs.size());
            vecrefs.reserve(oth.vecrefs.size());
            par2vec.reserve(oth.par2vec.size());
            u_vecrefs.reserve(oth.u_vecrefs.size());
            old2vals.reserve(oth.old2vals.size());
            if (oth.cav != nullptr) {
                cav = vecrefs.try_emplace(std::make_shared<Vector>(*(oth.cav)),
                        oth.vecrefs.at(oth.cav)).first->first;
                oth2this.try_emplace(oth.cav, cav);
            }
            for (const auto &p : oth.vecrefs)
                oth2this.try_emplace(p.first,
                        vecrefs.try_emplace(std::make_shared<Vector>(
                                *(p.first)),
                            p.second).first->first);
            for (const auto &p : oth.par2vec)
                par2vec.try_emplace(p.first, oth2this.at(p.second));
            for (const auto &p : oth.u_vecrefs)
                u_vecrefs.try_emplace(oth2this.at(p.first), p.second);
            for (const auto &p : oth.old2vals)
                old2vals.try_emplace(oth2this.at(p.first), p.second);
            return *this;
        }

        /** Move assignment operator. */
        ProjectionBase<T, F, V, Hash, EqualTo, LessV>& operator=(
                ProjectionBase<T, F, V, Hash, EqualTo, LessV> &&oth) noexcept {
            dims_ = oth.dims_;
            famsz_ = oth.famsz_;
            u_cps = oth.u_cps;
            u_idx = oth.u_idx;
            u_avg = oth.u_avg;
            u_pars = oth.u_pars;
            u_vecs = oth.u_vecs;
            cav = oth.cav;
            newvals = std::move(oth.newvals);
            par2vec = std::move(oth.par2vec);
            old2vals = std::move(oth.old2vals);
            u_vecrefs = std::move(oth.u_vecrefs);
            vecrefs = std::move(oth.vecrefs);
            return *this;
        }

        /** Begins an update.
         *
         * Prepares the Projection for an update, such as adding vectors,
         * mapping new parameters, and adding components to vectors.
         * Must be called before using update(). Once done with this update,
         * end_update() must be called, to actually perform the update.
         *
         * All updates involving specific component indices of the vectors in
         * the projection must either all use all the begin_update() - update()
         * - end_update() procedure or use all the update() (with four
         * parameters) function.
         *
         * @param idx  index of the first component to add.
         * @param cps  number of components to add.
         */
        virtual void begin_update(const std::size_t idx, const std::size_t cps) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (cps == 0)
                throw std::out_of_range(
                        "Wheelwright::Projection::begin_update: cps must be "
                        "positive.");
            if (u_cps != 0)
                throw std::out_of_range(
                        "Wheelwright::Projection::begin_update: begin_update() "
                        "called before the active update was ended with "
                        "end_update().");
#endif
            newvals.clear();
            old2vals.clear();
            old2vals.reserve(vecrefs.size());
            u_cps =  cps;
            u_idx = idx;
            u_avg = 0;
            u_pars = 0;
            u_vecrefs.clear();
            u_vecrefs.reserve(vecrefs.size());
            u_vecrefs.insert(vecrefs.begin(), vecrefs.end());
            if (! vecrefs.empty())
                u_vecs = vecrefs.size();
            else
                u_vecs = 1; // for the catch-all
        }

        /** Returns a const iterator to the first vector. */
        const_iterator cbegin() const {
            return const_iterator(vecrefs.cbegin());
        }

        /** Returns a const iterator past the last vector. */
        const_iterator cend() const { return const_iterator(vecrefs.cend()); }

        /** Empties the projection, resetting it as it was immediately after the
         * ctor completed.
         */
        virtual void clear() {
            dims_ = u_cps = u_idx = u_pars = u_vecs = 0;
            u_avg =  0;
            cav = nullptr;
            newvals.clear();
            par2vec.clear();
            old2vals.clear();
            u_vecrefs.clear();
            vecrefs.clear();
        }

        /** Returns the number of vectors that @a par is mapped to (0 or 1).
         *
         * @param par  the parameter to lookup.
         */
        std::size_t count(const T &par) const {
            return par2vec.count(par);
        }

        /** Returns the number of components of the vectors. */
        std::size_t dims() const { return dims_; }

        /** Sets the number of components of the vectors.
         *
         * @param s     the new number of components in the vectors.
         */
        void dims(const std::size_t s) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (vecrefs.empty())
                throw std::out_of_range("Wheelwright::Projection::dims: "
                        "setting the number of dimensions cannot be decreased");
            if (dims_ > s)
                throw std::out_of_range("Wheelwright::Projection::dims: the "
                         "number of dimensions cannot be decreased");
#endif
            dims_ = s;
        }

        /** Removes the catch-all vector from the projection.
         *
         * This function should be used when there is absolute certainty that no
         * additional parameters will be mapped to vectors in the projection,
         * i.e., when the number of mapped parameters equals the size of the
         * family of functions.
         */
        void drop_catchall() {
            vecrefs.erase(cav);
            cav = nullptr;
        }

        /** Ends an update.
         *
         * Performs the actual updating of the Projection.
         *
         * This function should be called after having completed a set of
         * updates to the projection using update().
         */
        virtual void end_update() {
#ifdef WHEELWRIGHT_INVARIANTS
            if (u_cps == 0)
                throw std::out_of_range(
                        "Wheelwright::Projection::end_update(): no update is "
                        "active.");
#endif
            // Check whether the Projection will no longer contain the catch-all
            // vector after the update, and if so, update the count of vectors,
            // so the avg is computed correctly.
            if (par2vec.size() + u_pars == famsz_)
                --u_vecs;
            u_avg /= u_vecs;
            if (const auto cdims {u_idx + u_cps}; cdims > dims_) {
                dims_ = cdims;
                // Start filling an empty Projection with the new vectors,
                // including, if needed, the catch-all vector.
                if (vecrefs.empty()) {
                    for (const auto &p : newvals) {
                        vec_p vec {vecrefs.try_emplace(
                                 built_from_scratch(p.first),
                                p.second.size()).first->first};
                        for (const T &par : p.second)
                            par2vec.emplace(par, vec);
                    }
                    if (famsz_ > u_pars) {
                        cav = built_from_scratch(0);
                        vecrefs.emplace(cav, famsz_ - par2vec.size());
                    }
                    // Reset the "guard" for active update.
                    u_cps = 0;
#ifdef WHEELWRIGHT_INVARIANTS
                    if (vecrefs.size() != u_vecs)
                        throw std::runtime_error(
                                "Wheelwright::Projection::end_update(): the "
                                "size of vecrefs is different than u_vecs."
                                );
                    if (par2vec.size() != u_pars)
                        throw std::runtime_error(
                                "Wheelwright::Projection::end_update(): the "
                                "size of par2vec is different than u_pars."
                                );
#endif
                    return;
                }
            }
            decltype(vecrefs) vecrefs_n;
            // Create and insert in the projection the new vectors that use the
            // existing ones as base (but not yet those that will use the
            // catch-all vector as base).
            for (const auto &p : old2vals) {
                std::size_t d {0};
                for (const auto &v : p.second) {
                    vec_p vec {vecrefs_n.try_emplace(
                            extended(*(p.first), v.first),
                            v.second.size()).first->first};
                    d += v.second.size();
                    for (const auto &par : v.second)
                        par2vec[par] = vec;
                }
                // The base vector may still have parameters mapped to it. In
                // this case, add the new components to it and insert it in
                // vecrefs_n.
                if (auto h = vecrefs.extract(p.first); h.mapped() > d) {
                    add_zeroes(h.key());
                    h.mapped() -= d;
                    vecrefs_n.insert(std::move(h));
                }
            }
            // Add new parameters mapped to new vectors that use the catch-all
            // vector as base, and then update the catch-all vector.
            if (cav != nullptr) {
                for (const auto &p : newvals) {
                    vec_p vec {vecrefs_n.try_emplace(
                            extended(*cav, p.first),
                            p.second.size()).first->first};
                    for (const T &par : p.second)
                        par2vec.emplace(par, vec);
                }
                if (auto cav_h {vecrefs.extract(cav)};
                        cav_h.mapped() > u_pars) {
                    add_zeroes(cav_h.key());
                    cav_h.mapped() -= u_pars;
                    // The use of std::move is from the example on
                    // http://en.cppreference.com/w/cpp/container/map/extract .
                    vecrefs_n.insert(std::move(cav_h));
                } else
                    cav = nullptr;
            }
            // The only vectors in vecrefs at this point are those that were not
            // split into new vectors. Add the new components to them and move
            // them to vecrefs_n.
            for (auto& p : vecrefs) {
                auto cp {p.first};
                add_zeroes(cp);
                vecrefs_n.emplace(cp, p.second);
            }
            vecrefs = std::move(vecrefs_n);
            // Reset the "guard" for active update.
            u_cps = 0;
        }

        /** Sets the upper bound to the family size.
         *
         * Throws an exception if the new upper bound is smaller than the
         * current number of parameters in the projection.
         *
         * Throws an exception if the current number of paraeters in the
         * projection is equal to the old upper bound to the family size.
         *
         * Discards the catch-all vector if the current number of parameters is
         * equal to the new upper bound.
         *
         * @param s     the new upper bound to the family size.
         */
        void famsz(const std::size_t s) {
            if (famsz_ == s)
                return;
            else if (famsz_ == par_size())
                throw std::runtime_error("Wheewlright::famsz: the current "
                        "number of vectors in the Projection is equal to the "
                        "old upper bound to the family size, thus the upper "
                        "bound cannot be changed.");
            else if (s < par_size())
                throw std::runtime_error("Wheewlright::famsz: the current "
                        "number of vectors in the Projection is larger than "
                        "the passed new upper bound.");
            else {
                if (par_size() == s)
                    drop_catchall();
                famsz_ = s;
            }
        }

        std::size_t famsz() const {
            return famsz_;
        }

        /** Returns an iterator to the first parameter mapped in the Projection.
         *
         * The order of the parameters is determined by @a Less.
         */
        auto par_cbegin() const {
            return boost::make_transform_iterator(par2vec.cbegin(),
                    MatteoUtils::KeyGetter<decltype(par2vec.cbegin())>());
        }

        /** Returns an iterator past the last parameter mapped in the
         * Projection.
         *
         * The order of the parameters is determined by @a Less.
         */
        auto par_cend() const {
            return boost::make_transform_iterator(par2vec.cend(),
                    MatteoUtils::KeyGetter<decltype(par2vec.cend())>());
        }

        /** Returns the number of parameters that are mapped to the vectors. */
        std::size_t par_size() const { return par2vec.size(); }

        /** Returns a const reference to the pointer to the vector that @a par
         * maps to.
         *
         * The function makes no range checking on @a par.
         *
         * @param par  the T parameter mapped to the vector to return.
         */
        vec_p const& operator[](const T par) const {
            if (auto it {par2vec.find(par)}; it == par2vec.end())
                return cav;
            else
                return it->second;
        }

        /** Returns the number of parameters mapped to this vector.
         *
         * The function makes no range checking on @a v.
         *
         * @param v  the vec_p whose number of references to return.
         */
         std::size_t refs(const vec_p v) { return vecrefs[v]; }

        /** Returns the number of vectors. */
        std::size_t size() const { return vecrefs.size(); }

        /** Updates the vector to which @a par is mapped to, or maps @a par to a
         * (potentially new) vector.
         *
         * This function must be called after begin_update() has been called.
         *
         * This function does not really perform the updates, rather it
         * registers it for future execution, which will be carried on by
         * end_update().
         *
         * @param par  the parameter to update.
         * @param val  the value for the components to add.
         **/
        void update(const T &par, V val) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (u_cps == 0)
                throw std::out_of_range(
                        "Wheelwright::Projection::update: update() called but "
                        "no update is active.");
#endif
            if (auto par_it {par2vec.find(par)}; par_it == par2vec.end()) {
                // "par" is a new parameter never seen before.
#ifdef WHEELWRIGHT_INVARIANTS
                // TODO: Wheelwright::impl::ProjectionBase::update: Why the +1?
                if (par2vec.size() + u_pars == famsz_ + 1)
                    throw std::runtime_error(
                            "Wheelwright::Projection::update: the number of "
                            "parameters mapped to vectors (i.e., the size of the "
                            "family of functions), has already reached the value "
                            "passed in the constructor.");
#endif
                ++u_pars;
                if (auto r {newvals.try_emplace(val,
                            typename decltype(newvals)::mapped_type {par})};
                        r.second == true) {
                    u_avg += val;
                    ++u_vecs;
                } else
                    r.first->second.push_back(par);
            } else {
                // par is in par2vec, we saw it before.

                // The expression used to construct r in the following test
                // first adds (if it doesn't exist) an entry in old2vals from
                // the Vector that par is mapped to to an empty map, and then
                // adds (if it doesn't exist) an entry in this map (which may
                // not be empty) from val to an empty std::vector.
                if (auto r {old2vals.try_emplace(par_it->second,
                        typename decltype(old2vals)::mapped_type
                        {}).first->second.try_emplace(val, std::vector<T> {})};
                        r.second == true) {
                    r.first->second.push_back(par);
                    u_avg += val;
                    ++u_vecs;
                } else
                    // TODO: Wheelwright::ProjectionBase::update: devise a test
                    // to stress this piece of code.
                    r.first->second.push_back(par);
                // Decrease the number of references that the old vector will
                // have after the update is complete, and if it reaches zero,
                // just decrease its count.
                if (auto u_it {u_vecrefs.find(par_it->second)};
                        u_it->second > 1)
                    u_it->second -= 1;
                else {
                    --u_vecs;
                    u_vecrefs.erase(u_it);
                }
            }
        }
}; // class ProjectionBase


template<typename T> class DiscreteBase {

    public:

        DiscreteBase(const DiscreteBase<T> &) = default;
        DiscreteBase(DiscreteBase<T> &&) = default;
        DiscreteBase<T>& operator=(const DiscreteBase<T> &) = default;
        DiscreteBase<T>& operator=(DiscreteBase<T> &&) = default;

    protected:
        template<typename Iterator> DiscreteBase(Iterator b,
                Iterator e, typename std::enable_if_t<!
                std::is_arithmetic<Iterator>::value>* = 0) : codomain(b, e) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (codomain.size() == 0)
                throw std::runtime_error("Wheelwright::DiscreteBase: "
                        "codomain cannot be empty.");
#endif
        }

        template<typename... J> DiscreteBase(J ...args,
                std::enable_if_t<sizeof...(J) != 0> = 0) : codomain {args...} {
#ifdef WHEELWRIGHT_INVARIANTS
            if (codomain.size() == 0)
                throw std::runtime_error("Wheelwright::DiscreteBase: "
                        "codomain cannot be empty.");
#endif
        }

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        void do_merge(DiscreteBase &oth) {
#ifdef WHEELWRIGHT_INVARIANTS
            if (oth.codomain.size() != codomain.size())
                throw std::runtime_error("Wheelwright::DiscreteBase::do_merge: "
                        "the two Projections do not have the same co-domain.");
            else
                for (auto vi {codomain.cbegin()}, oi {oth.codomain.cbegin()};
                        vi != codomain.cend(); ++vi, ++oi) {
                    if (*vi != *oi)
                        throw std::runtime_error(
                                "Wheelwright::DiscreteBase::do_merge: the two "
                                "Projections do not have the same co-domain.");
                }
#endif
            std::ignore = oth;
        }

        std::set<T> codomain;
};


class OneDrawBase {

    public:
        /** Set the seed of the generator of Rademacher variables to @a s.
         *
         * This function is used in the tests to make them replicable.
         */
        void seed(std::size_t s) {
            gen.seed(s);
        }

    protected:
        OneDrawBase() : gen {std::random_device{}()} {}

        void sample_rademacher(const std::size_t cps) {
            u_rade.reserve(cps);
            for (std::size_t i {0}; i < cps; ++i)
                u_rade.push_back((dis(gen)) ? Rademacher::POS :
                        Rademacher::NEG);
        }

        // We can rely on the implicitly-generated copy/move constructors and
        // assignment operators.

        std::vector<Rademacher> u_rade;

    private:
        std::mt19937 gen;
        std::uniform_int_distribution<signed char> dis {0, 1};
}; // class OneDrawBase

} // namespace Wheelwright::impl

#endif // WHEELWRIGHT_PROJECTION_IMPL_H
