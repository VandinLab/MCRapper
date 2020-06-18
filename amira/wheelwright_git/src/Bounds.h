/** @file Bounds.h
 * Upper bounds to the Empirical Rademacher Average (ERA) and to the maximum
 * deviation.
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

#ifndef WHEELWRIGHT_BOUNDS_H_
#define WHEELWRIGHT_BOUNDS_H_

#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "matteoutils/Iterators.h"
#include "matteoutils/NumericTools.h"

#include "Bounds_impl.h"

/** General namespace for the library.
 *
 * All functions and classes belong to this namespace.
 *
 * Implementation details are in the non-documented  ``Wheelwright::impl``
 * namespace.
 */
namespace Wheelwright {

/** Specifies whether to the compute upper bound is for the maximum deviation or
 * the maximum absolute deviation of sample averages from their expectations.
 */
enum class MaxDevType: unsigned char {
    /** "Plain" maximum deviation. */
    PLAIN = 1,
    /** Absolute maximum deviation. */
    ABSOL = 2
};

/** Method to bound the Rademacher average.
 */
enum class RaUBMethod : unsigned char {
    /** Use the Empirical Rademacher Average and McDiarmid's bounded difference
     * inequality.
     */
    MCDIARMID = 1,
    /** Use the Empirical Rademacher Average and the tail inequality for
     * self-bounding functions.
     */
    SELFBOUNDING = 2,
    /** Use the one-draw approach. */
    ONEDRAW = 3
};

/** Returns a string representation of @a d. */
inline std::string to_string_RaUBMethod(const RaUBMethod d) {
    switch (d) {
        case RaUBMethod::MCDIARMID:
            return "MCDIARMID";
        case RaUBMethod::SELFBOUNDING:
            return "SELFBOUNDING";
        case RaUBMethod::ONEDRAW:
            return "ONEDRAW";
        default:
            throw std::out_of_range("Wheelwright::to_string_RaUBMethod: "
                    "Reached code supposed to be unreachable");
    }
}


/** Method to compute the upper bound to the Empirical Rademacher Average (ERA).
 */
enum class EraUBMethod : unsigned char {
    /** Use Anguita et al's bound. See era_bound_AGOR(). */
    AGOR = 1,
    /** Use Anguita et al's bound specialized for vectors with components from
     * a finite discrete domain. See era_bound_AGOR_FD(). */
    AGOR_FD = 2,
    /** Use Massart's Lemma. See era_bound_MASSART(). */
    MASSART = 3,
    /** Use Riondato & Upfal's bound. See era_bound_RU(). */
    RU = 4,
    /** Use Riondato & Upfal's bound specialized for vectors with components from
     * a finite discrete domain. See era_bound_RU_FD(). */
    RU_FD =5
};

/** Returns a string representation of the EraUBMethod @a e. */
inline std::string to_string_EraUBMethod(const EraUBMethod e) {
    switch (e) {
        case EraUBMethod::AGOR:
            return "AGOR";
        case EraUBMethod::AGOR_FD:
            return "AGOR_FD";
        case EraUBMethod::MASSART:
            return "MASSART";
        case EraUBMethod::RU:
            return "RU";
        case EraUBMethod::RU_FD:
            return "RU_FD";
        default:
            throw std::out_of_range("Wheelwright::to_string_EraUBMethod: "
                    "Reached code supposed to be unreachable");
    }
}

/** Returns Anguita et al's ERA upper bound (random-access-iterator variant).
 *
 * Returns an upper bound to the empirical Rademacher average (ERA) computed
 * using a generalized version of Anguita et al.'s method presented in Sect. III
 * of D. Anguita, A. Ghio, L. Oneto, and S. Ridella, "A Deep Connection Between
 * the Vapnik–Chervonenkis Entropy and the Rademacher Complexity", IEEE Trans.
 * Neural Netw. Learn. Sys., Vol. 25, No. 12, Dec. 2014, pp. 2202--2211.
 *
 * See also the era_bound_RU function, which sometimes offers a tighter bound.
 *
 * This function is available only when the `WHEELWRIGHT_NLOPT` preprocessor
 * directive is defined.
 *
 * @tparam  Iterator : an Iterator type
 *
 * @param begin   beginning of a sequence of ProjVec's
 * @param end     past-end of a sequence of ProjVec's
 * @param dim     the dimension (no. of components) of the vectors
 * @param ftol    the tolerance for early stopping of the optimization. Default:
 *              1e-7.
 * @param unused  enabler for the OpenMP version when @a Iterator is a random
 *                access iterator.
 */
template<class Iterator> double era_bound_AGOR(Iterator begin, Iterator end,
        const std::size_t dim, const double ftol = 1e-7,
        std::enable_if_t<
            MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
#ifndef WHEELWRIGHT_NLOPT
#error "The optimization-based bounds to the empirical Rademacher average are not available because Wheelwright is being compiled without NLopt support ('WHEELWRIGHT_NLOPT' is not defined)."
    std::ignore = begin;
    std::ignore = end;
    std::ignore = dim;
    std::ignore = ftol;
#endif // WHEELWRIGHT_NLOPT
    const auto dist {std::distance(begin,end)};
    std::vector<double> objSumExponents(dist, 0);
    #pragma omp parallel for
    for (auto i = decltype(dist){0}; i < dist; ++i)
        objSumExponents[i] = impl::elltwosquared(*(begin + i)) / 2;
    return impl::optimize(impl::objective_AGORRU, &objSumExponents, ftol) / dim;
}

/** Returns Anguita et al's ERA upper bound (non-random-access-iterator
 * variant).
 *
 * See the documentation for the random-access-iterator variant.
 */
template<class Iterator> double era_bound_AGOR(Iterator begin, Iterator end,
        const std::size_t dim, const double ftol = 1e-7,
        std::enable_if_t<
            ! MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
#ifndef WHEELWRIGHT_NLOPT
#error "The optimization-based bounds to the empirical Rademacher average are not available because Wheelwright is being compiled without NLopt support ('WHEELWRIGHT_NLOPT' is not defined)."
    std::ignore = begin;
    std::ignore = end;
    std::ignore = dim;
    std::ignore = ftol;
#endif // WHEELWRIGHT_NLOPT
    std::vector<double> objSumExponents;
    for (auto i = begin; i != end; ++i)
        objSumExponents.emplace_back(impl::elltwosquared(*i) / 2);
    return impl::optimize(impl::objective_AGORRU, &objSumExponents, ftol) / dim;
}

/** Returns Anguita et al's ERA upper bound for family of functions with a
 * finite, discrete co-domain. (random-access-iterator variant)
 *
 * Returns an upper bound to the empirical Rademacher average (ERA) computed
 * using Anguita et al.'s bound, specialized for functions with a finite,
 * discrete co-domain.
 *
 * This function is available only when the `WHEELWRIGHT_NLOPT` preprocessor
 * directive is defined.
 *
 * @tparam VecIter the type for the sequence of vectors.
 * @tparam DomIter the type for the sequence of vectors.
 *
 * @param [begin,end)  the sequence of ProjVec's.
 * @param dim          the dimension (no. of components) of the vectors.
 */
template<typename VecIter, typename DomIter> double era_bound_AGOR_FD(
        VecIter begin, VecIter end, const std::size_t dim, DomIter dom_begin,
        DomIter dom_end, const double ftol = 1e-7, std::enable_if_t<
        MatteoUtils::is_random_access_iterator<VecIter>::value>* = 0) {
    // The arguments for the cosh functions. In this case (AGOR), just the
    // elements of the domain.
    // We assume that the type of the entries of the vectors can be
    // automatically cast to double.
    std::vector<double> args {dom_begin, dom_end};
    const auto dist {std::distance(begin, end)};
    std::vector<std::vector<std::size_t>> counts(dist);
    #pragma omp parallel for
    for (auto i = decltype(dist){0}; i < dist; ++i) {
        counts[i] = {(begin + i)->begin(), (begin + i)->end()};
#ifdef WHEELWRIGHT_INVARIANTS
        if (counts[i].size() != args.size()) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_AGOR_FD: found a vector with "
                << counts[i].size() << " counts but the size of the domain "
                "is "  << args.size() << ".";
            throw std::runtime_error(ss.str());
        }
        std::size_t c {0};
        #pragma omp parallel for reduction(+: c)
        for (std::size_t j = 0; j < counts[i].size(); ++j) {
            c += counts[i][j];
        }
        if (c > dim) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_AGOR_FD: the sum of the counts for a "
                "vector is " << c << " but it cannot not be larger than " << dim
                << "the passed number of components.";
            throw std::runtime_error(ss.str());
        }
#endif
    }
    auto args_and_counts {std::pair(args, counts)};
    return impl::optimize(impl::objective_AGORRU_FD, &args_and_counts, ftol) /
        dim;
}

/** Returns Anguita et al's ERA upper bound for family of functions with a
 * finite, discrete co-domain (non-random-access-iterator variant).
 *
 * See the documentation for the random-access-iterator variant.
 */
template<typename VecIter, typename DomIter> double era_bound_AGOR_FD(
        VecIter begin, VecIter end, const std::size_t dim, DomIter dom_begin,
        DomIter dom_end, const double ftol = 1e-7, std::enable_if_t<
        ! MatteoUtils::is_random_access_iterator<VecIter>::value>* = 0) {
    // The arguments for the cosh functions. In this case (AGOR), just the
    // elements of the domain.
    // We assume that the type of the entries of the vectors can be
    // automatically cast to double.
    std::vector<double> args {dom_begin, dom_end};
    std::vector<std::vector<std::size_t>> counts;
    for (; begin != end; ++begin) {
        counts.emplace_back(begin->begin(), begin->end());
#ifdef WHEELWRIGHT_INVARIANTS
        if (counts[counts.size() - 1].size() != args.size()) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_AGOR_FD: found a vector with "
                << counts[counts.size() - 1].size() << " counts but the size "
                "of the domain is "  << args.size() << ".";
            throw std::runtime_error(ss.str());
        }
        std::size_t c {0};
        #pragma omp parallel for reduction(+: c)
        for (std::size_t j = 0; j < counts[counts.size() -1].size(); ++j) {
            c += counts[counts.size() - 1][j];
        }
        if (c > dim) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_AGOR_FD: the sum of the counts for a "
                "vector is " << c << " but it cannot not be larger than " << dim
                << "the passed number of components.";
            throw std::runtime_error(ss.str());
        }
#endif
    }
    auto args_and_counts {std::pair(args, counts)};
    return impl::optimize(impl::objective_AGORRU_FD, &args_and_counts, ftol) /
        dim;
}

/** Returns the ERA upper bound from Massart's Lemma.
 *
 * Returns an upper bound to the empirical Rademacher average (ERA) computed
 * using Massart's Lemma on the passed values describing a projection.
 *
 * @tparam T  the numeric type of @a maxelltwosquared.
 *
 * @param maxelltwo  the maximum @f$\ell_2@f$-norm of a vector in the projection.
 * @param vecs       the number of vectors in the projection.
 * @param dim       the dimension (no. of components) of the vectors.
 */
template<typename T> double era_bound_MASSART(const T maxelltwo,
        const std::size_t vecs, const std::size_t dim,
        std::enable_if_t<std::is_arithmetic<T>::value>* = 0) {
    return maxelltwo * std::sqrt(2 * std::log(vecs)) / dim;
}

/** Returns the ERA upper bound from Massart's Lemma (random-access-iterator
 * variant).
 *
 * Returns an upper bound to the empirical Rademacher average (ERA) computed
 * using Massart's Lemma on the projection containing the dimensional vectors in
 * [@a begin, @a end), each with @a size components. A reference for Massart's
 * Lemma is Lemma 26.8 in S. Shalev-Shwartz and S. Ben-David, "Understanding
 * Machine Learning", Cambridge University Press, 2014.
 *
 * @tparam Iterator  an Iterator type.
 *
 * @param begin  beginning of a sequence of ProjVec's.
 * @param end    end of a sequence of ProjVecs.
 * @param size   the sample size
 */
template<typename Iterator> double era_bound_MASSART(Iterator begin,
        Iterator end, const std::size_t size, std::enable_if_t<
            MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
    const auto vecs {std::distance(begin,end)};
    return era_bound_MASSART(std::sqrt(impl::maxelltwosquared(begin, end)),
            vecs, size);
}

/** Returns the ERA upper bound from Massart's Lemma (non-random-access-iterator
 * variant).
 *
 * See the documentation for the random-access-iterator variant.
 */
template<typename Iterator> double era_bound_MASSART(Iterator begin,
        Iterator end, const std::size_t size,
        std::enable_if_t<MatteoUtils::is_iterator<Iterator>::value>* = 0,
        std::enable_if_t<
            ! MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
    Iterator bc {begin};
    Iterator endc {end};
    return era_bound_MASSART(std::sqrt(impl::maxelltwosquared(begin, end)),
            std::distance(bc, endc), size);
}

/** Returns Riondato and Upfal's upper bound to the ERA (random-access-iterator
 * variant).
 *
 * Returns an upper bound to the empirical Rademacher average (ERA) computed
 * according to the formula in Thm. 3 of Riondato M. and Upfal E., "Mining
 * Frequent Itemsets through Progressive Sampling with Rademacher Averages", ACM
 * KDD'15 (see the extended version at
 * http://matteo.rionda.to/papers/RiondatoUpfal-FrequentItemsetsSamplingRademacher-KDD.pdf
 * .)
 *
 * The bound returned by this function is often much tighter than the one
 * returned by era_bound_AGOR, especially for functions taking values greater
 * than one.
 *
 * This function is available only when the `WHEELWRIGHT_NLOPT` preprocessor
 * directive is defined.
 *
 * @tparam Iterator  an Iterator type.
 *
 * @param begin   beginning of a sequence of ProjVec's.
 * @param end     end of a sequence of ProjVecs.
 * @param dim     the dimension (no. of components) of the vectors.
 * @param ftol    the tolerance for early stopping of the optimization. Default:
 *              1e-7.
 * @param unused  enabler for the OpenMP version when @a Iterator is a random
 *              access iterator.
 */
template<class Iterator> double era_bound_RU(Iterator begin, Iterator end,
        const std::size_t dim, const double ftol = 1e-7,
        std::enable_if_t<
            MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
#ifndef WHEELWRIGHT_NLOPT
#error "The optimization-based bounds to the empirical Rademacher average are not available because Wheelwright is being compiled without NLopt support ('WHEELWRIGHT_NLOPT' is not defined)."
    std::ignore = begin;
    std::ignore = end;
    std::ignore = dim;
    std::ignore = ftol;
#endif // WHEELWRIGHT_NLOPT
    const double exponentDenom {2 * std::pow(dim, 2)};
    const auto dist {std::distance(begin, end)};
    std::vector<double> objSumExponents(dist, 0);
    #pragma omp parallel for
    for (auto i = decltype(dist){0}; i < dist; ++i)
        objSumExponents[i] = impl::elltwosquared(*(begin + i)) / exponentDenom;
    return impl::optimize(impl::objective_AGORRU, &objSumExponents, ftol);
}

/** Returns Riondato and Upfal's upper bound to the ERA
 * (non-random-access-iterator variant).
 *
 * See the documentation for the random-access-iterator variant.
 */
template<class Iterator> double era_bound_RU(Iterator begin, Iterator end,
        const std::size_t dim, const double ftol = 1e-7,
        std::enable_if_t<
            ! MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
#ifndef WHEELWRIGHT_NLOPT
#error "The optimization-based bounds to the empirical Rademacher average are not available because Wheelwright is being compiled without NLopt support ('WHEELWRIGHT_NLOPT' is not defined)."
    std::ignore = begin;
    std::ignore = end;
    std::ignore = dim;
    std::ignore = ftol;
#endif // WHEELWRIGHT_NLOPT
    const double exponentDenom {2 * std::pow(dim, 2)};
    std::vector<double> objSumExponents;
    for (; begin != end; ++begin)
        objSumExponents.push_back(impl::elltwosquared(*begin) / exponentDenom);
    return impl::optimize(impl::objective_AGORRU, &objSumExponents, ftol);
}

/** Returns Riondato and Upfal's ERA upper bound for family of functions with a
 * finite, discrete co-domain (random-access-iterator variant).
 *
 *
 * Returns an upper bound to the empirical Rademacher average (ERA) computed
 * using Riondato and Upfal's bound, specialized for functions with a finite,
 * discrete co-domain.
 *
 * This function is available only when the `WHEELWRIGHT_NLOPT` preprocessor
 * directive is defined.
 *
 * @tparam VecIter the type for the sequence of vectors.
 * @tparam DomIter the type for the sequence of vectors.
 *
 * @param [begin,end)  the sequence of ProjVec's.
 * @param dim          the dimension (no. of components) of the vectors.
 */
template<typename VecIter, typename DomIter> double era_bound_RU_FD(
        VecIter begin, VecIter end, const std::size_t dim, DomIter dom_begin,
        DomIter dom_end, const double ftol = 1e-7, std::enable_if_t<
        MatteoUtils::is_random_access_iterator<VecIter>::value>* = 0) {
    // The arguments for the cosh functions. In this case (RU), the
    // elements of the domain divided by the number of dimensions.
    std::vector<double> args;
    for (; dom_begin != dom_end; ++dom_begin)
        args.emplace_back(((double) *dom_begin) / dim);
    const auto dist {std::distance(begin, end)};
    std::vector<std::vector<std::size_t>> counts(dist);
    #pragma omp parallel for
    for (auto i = decltype(dist){0}; i < dist; ++i) {
        counts[i] = {(begin + i)->begin(), (begin + i)->end()};
#ifdef WHEELWRIGHT_INVARIANTS
        if (counts[i].size() != args.size()) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_RU_FD: found a vector with "
                << counts[i].size() << " counts but the size of the domain "
                "is "  << args.size() << ".";
            throw std::runtime_error(ss.str());
        }
        std::size_t c {0};
        #pragma omp parallel for reduction(+: c)
        for (std::size_t j = 0; j < counts[i].size(); ++j) {
            c += counts[i][j];
        }
        if (c > dim) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_RU_FD: the sum of the counts for a "
                "vector is " << c << " but it cannot not be larger than " << dim
                << "the passed number of components.";
            throw std::runtime_error(ss.str());
        }
#endif
    }
    auto args_and_counts {std::pair(args, counts)};
    return impl::optimize(impl::objective_AGORRU_FD, &args_and_counts, ftol);
}

/** Returns Riondato and Upfal's ERA upper bound for family of functions with a
 * finite, discrete co-domain (non-random-access-iterator variant).
 *
 * See the documentation for the random-access-iterator variant.
 */
template<typename VecIter, typename DomIter> double era_bound_RU_FD(
        VecIter begin, VecIter end, const std::size_t dim, DomIter dom_begin,
        DomIter dom_end, const double ftol = 1e-7, std::enable_if_t<
        ! MatteoUtils::is_random_access_iterator<VecIter>::value>* = 0) {
    // The arguments for the cosh functions. In this case (RU), the
    // elements of the domain divided by the number of dimensions.
    std::vector<double> args;
    for (; dom_begin != dom_end; ++dom_begin)
        args.emplace_back(((double) *dom_begin) / dim);
    std::vector<std::vector<std::size_t>> counts;
    for (; begin != end; ++begin) {
        counts.emplace_back(begin->begin(), begin->end());
#ifdef WHEELWRIGHT_INVARIANTS
        if (counts[counts.size() - 1].size() != args.size()) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_RU_FD: found a vector with "
                << counts[counts.size() - 1].size() << " counts but the size "
                "of the domain is "  << args.size() << ".";
            throw std::runtime_error(ss.str());
        }
        std::size_t c {0};
        #pragma omp parallel for reduction(+: c)
        for (std::size_t j = 0; j < counts[counts.size() - 1].size(); ++j) {
            c += counts[counts.size() - 1][j];
        }
        if (c > dim) {
            std::stringstream ss;
            ss <<"Wheelwright::era_bound_RU_FD: the sum of the counts for a "
                "vector is " << c << " but it cannot not be larger than " << dim
                << "the passed number of components.";
            throw std::runtime_error(ss.str());
        }
#endif
    }
    auto args_and_counts {std::pair(args, counts)};
    return impl::optimize(impl::objective_AGORRU_FD, &args_and_counts, ftol);
}

/** Returns an upper bound to the Empirical Rademacher Average computed using
 * the specified method.
 *
 * Variant using iterators.
 *
 * @tparam Iterator  an Iterator type.
 * @tparam Args      variadic template type for additional arguments.
 *
 * @param begin  beginning of a sequence of ProjVec's.
 * @param end    end of a sequence of ProjVecs.
 * @param size   the sample size.
 * @param e      the EraUBMethod method to use for computing the returned value.
 * @param args   additional arguments.
 */
template<typename Iterator, typename... Args> double era_bound(Iterator begin,
        Iterator end, const std::size_t size, const EraUBMethod e,
        Args... args) {
    // This function must be defined after the era_bound_* ones or the code
    // won't compile.
    if constexpr (sizeof...(args) < 2) {
        switch (e) {
            case EraUBMethod::AGOR:
                return era_bound_AGOR(begin, end, size, args...);
                break;
            case EraUBMethod::MASSART:
                return era_bound_MASSART(begin, end, size);
                break;
            case EraUBMethod::RU:
                return era_bound_RU(begin, end, size, args...);
                break;
            default:
                throw std::out_of_range("Wheelwright::era_bound: e must be one "
                        "of AGOR, MASSART, or RU.");
        }
    } else {
        switch (e) {
            case EraUBMethod::AGOR_FD:
                return era_bound_AGOR_FD(begin, end, size, args...);
                break;
            case EraUBMethod::RU_FD:
                return era_bound_RU_FD(begin, end, size, args...);
                break;
            default:
                std::stringstream ss;
                ss << "Wheelwright::era_bound: the passed ERAUBMethod '"
                    << to_string_EraUBMethod(e)
                    << "' is not one of AGOR_FD or RU_FD.";
                throw std::out_of_range(ss.str());
        }
    }
}

/** Returns an upper bound to the Empirical Rademacher Average computed using
 * the specified method.
 *
 * Variant using Projection.
 *
 * @tparam P      the type of the @a proj argument, i.e., a specialization of
 *                Projection.
 * @tparam Args   variadic template type for additional arguments.
 *
 * @param proj   the Projection containing the Vectors for the functions on the
 *               sample.
 * @param e      the EraUBMethod method to use for computing the returned value.
 */
template<typename P, typename... Args> double era_bound(const P &proj,
        const EraUBMethod e, Args... args) {
    //if constexpr(impl::has_domain_begin<
            //typename std::iterator_traits<
            //typename P::const_iterator>::value_type>::value) {
        if (e == EraUBMethod::AGOR_FD || e == EraUBMethod::RU_FD)
            return era_bound(proj.cbegin(), proj.cend(), proj.dims(), e,
                    proj.cbegin()->domain_begin(), proj.cbegin()->domain_end(),
                    args...);
        else
            return era_bound(proj.cbegin(), proj.cend(), proj.dims(), e, args...);
    //} else
        //return era_bound(proj.cbegin(), proj.cend(), proj.dims(), e, args...);
}

/** Returns an upper bound to the Rademacher average.
 *
 * @param c     if @a r is Wheelwright::RaUBMethod::ONEDRAW, the ODARA,
 *              otherwise, the empirical Rademacher average (or an upper bound).
 * @param size  the sample size.
 * @param eta   the probability of error (i.e., 1 - confidence).
 * @param span  the width of the co-domain of the functions (or an upper bound).
 * @param r     the RaUBMethod method to use for computing the returned value.
 */
inline double ra_bound(const double c, const std::size_t size,
        const double eta, const double span, const RaUBMethod r) {
    if (r == Wheelwright::RaUBMethod::SELFBOUNDING)
        return impl::ra_bound_SELFBOUNDING(c, size, eta, span);
    else
        return impl::ra_bound_MCDIARMIDONEDRAW(c, size, eta, span);
}

/** Returns an upper bound to the Rademacher average.
 *
 * Variant using Projection.
 *
 * @tparam P  the type of the @a proj argument, i.e., a specialization of
 *            Projection.
 *
 * @param proj  the Projection containing the Vectors for the functions on the
 *              sample.
 * @param eta   the probability of error (i.e., 1 - confidence).
 * @param span  the width of the co-domain of the functions (or an upper bound).
 * @param r     the RaUBMethod method to use for computing the returned value.
 * @param e     the EraUBMethod method to use for computing the upper bound to
 *              the empirical Rademacher average. Ignored if @a r is
 *              Wheelwright::RaUBMethod::ONEDRAW.
 */
template<typename P> double ra_bound(const P &proj, const double eta,
        const double span, const RaUBMethod r, const EraUBMethod e) {
    return ra_bound(proj.cbegin(), proj.cend(), proj.dims(), eta, span, r, e);
}

/** Returns an upper bound to the Rademacher average.
 *
 * Variant using iterators.
 *
 * @tparam Iterator  the type of the @a begin and @a end arguments, behaving
 *                   like an iterator.
 *
 * @param [begin,end) the sequence of Vectors or numeric values representing the
 *                    appropriate quanties used to compute the upperbound to the
 *                    Rademacher average.
 * @param eta   the probability of error (i.e., 1 - confidence).
 * @param span  the width of the co-domain of the functions (or an upper bound).
 * @param r     the RaUBMethod method to use for computing the returned
 *              value.
 * @param e     the EraUBMethod method to use for computing the upper bound to
 *              the empirical Rademacher average. Ignored if @a r is
 *              Wheelwright::RaUBMethod::ONEDRAW.
 */
template<typename Iterator> double ra_bound(const Iterator begin,
        const Iterator end, const std::size_t size, const double eta,
        const double span, const RaUBMethod r, const EraUBMethod e) {
    if constexpr
        (impl::has_onedrawsum<typename
         std::iterator_traits<Iterator>::value_type>::value)
        if (r ==  Wheelwright::RaUBMethod::ONEDRAW)
            return impl::ra_bound_MCDIARMIDONEDRAW(odara(begin, end, size),
                    size, eta, span);
    return ra_bound(era_bound(begin, end, size, e), size, eta, span, r);
}

/** Returns an upper bound to the maximum deviation of sample averages from
 * their expectation for a family of functions.
 *
 * Variant using a numeric value.
 *
 * @param ra    the Rademacher average (or an upper bound to it).
 * @param size  the sample size.
 * @param eta   the probability of error (i.e., 1 - confidence).
 * @param span  the width of the co-domain of the functions (or an upper bound).
 * @param t     the kind of tail bound that is requested.
 */
inline double dev_bound(const double ra, const std::size_t size,
        const double eta, const double span, const MaxDevType t) {
    const double c {(t == MaxDevType::PLAIN) ? 0 : std::log(2)};
    return 2 * ra + span * std::sqrt((c - std::log(eta)) / (2 * size));
}

/** Returns an upper bound to the maximum deviation of sample averages from
 * their expectation for a family of functions.
 *
 * Variant using Projection.
 *
 * @tparam P  the type of the @a proj argument, i.e., a specialization of
 *            Projection.
 *
 * @param proj  the Projection containing the Vectors for the functions on the
 *              sample.
 * @param eta   the probability of error (i.e., 1 - confidence).
 * @param span  the width of the co-domain of the functions (or an upper bound).
 * @param r     the RaUBMethod method to use for computing the returned value.
 * @param e     the EraUBMethod method to use for computing the upper bound to
 *              the empirical Rademacher average.
 * @param t     the kind of tail bound that is requested.
 */
template<typename P> double dev_bound(const P &proj, const double eta,
        const double span, const RaUBMethod r, const EraUBMethod e,
        const MaxDevType t) {
    return dev_bound(proj.cbegin(), proj.cend(), proj.dims(), eta, span, r, e,
            t);
}

/** Returns an upper bound to the maximum deviation of sample averages from
 * their expectation for a family of functions.
 *
 * Variant using iterators.
 *
 * @tparam Iterator  the type of the @a begin and @a end arguments, behaving
 *                   like an iterator.
 *
 * @param [begin,end) the sequence of Vectors or numeric values representing the
 *                    appropriate quanties used to compute the upperbound to the
 *                    Rademacher average.
 * @param size        the sample size.
 * @param eta         the probability of error (i.e., 1 - confidence).
 * @param span        the width of the co-domain of the functions (or an upper
 *                    bound).
 * @param r           the RaUBMethod method to use for computing the returned
 *                    value.
 * @param e           the EraUBMethod method to use for computing the upper
 *                    bound to the empirical Rademacher average.
 * @param t           the kind of tail bound that is requested.
 */
template<typename Iterator> double dev_bound(const Iterator begin,
        const Iterator end, const double size, const double eta,
        const double span, const RaUBMethod r, const EraUBMethod e,
        const MaxDevType t) {
    const double adj_eta {(r == RaUBMethod::SELFBOUNDING || t ==
            MaxDevType::ABSOL) ?  eta / 2 : eta};
    const double ra {ra_bound(begin, end, size, adj_eta, span, r, e)};
    return dev_bound(ra, size,
            (r == Wheelwright::RaUBMethod::SELFBOUNDING) ? adj_eta: eta, span,
            t);
}


/** Compute the ODARA (random-access-iterator variant).
 *
 * @tparam Iterator  the type of the @a begin and @a end arguments, behaving
 *                   like an iterator.
 *
 * @param [begin,end) the sequence of Vectors or numeric values representing the
 *                    appropriate quanties used to compute the upperbound to the
 *                    Rademacher average.
 * @param size  the sample size.
 */
template<typename Iterator> double odara(const Iterator begin,
        const Iterator end, const std::size_t size, std::enable_if_t<
            MatteoUtils::is_random_access_iterator<Iterator>::value>*) {
    auto min {std::numeric_limits<double>::max()};
    auto max {std::numeric_limits<double>::lowest()};
    #pragma omp parallel for reduction(max : max) reduction(min: min)
    for (decltype(std::distance(begin,end)) i = 0;
            i < std::distance(begin, end); ++i) {
        const auto s {onedrawsum(*(begin + i))};
        if (min > s)
            min = s;
        if (max < s)
            max = s;
    }
    return (max - min) / (2 * size);
}

/** Compute the ODARA (non-random-access-iterator variant).
 *
 * @tparam Iterator  the type of the @a begin and @a end arguments, behaving
 *                   like an iterator.
 *
 * @param [begin,end) the sequence of Vectors or numeric values representing the
 *                    appropriate quanties used to compute the upperbound to the
 *                    Rademacher average.
 * @param size  the sample size.
 */
template<typename Iterator> double odara(const Iterator begin,
        const Iterator end, const std::size_t size,
        std::enable_if_t< !
            MatteoUtils::is_random_access_iterator<Iterator>::value>*) {
    auto min {std::numeric_limits<double>::max()};
    auto max {std::numeric_limits<double>::lowest()};
    for (auto it {begin}; it != end; ++it) {
        const auto s {onedrawsum(*it)};
        if (min > s)
            min = s;
        if (max < s)
            max = s;
    }
    return (max - min) / (2 * size);
}

/** Compute the ODARA (Projection variant).
 *
 * @tparam P  the type of the @a proj argument, i.e., a specialization of
 *            Projection.
 *
 * @param proj   the Projection containing the Vectors for the functions on the
 *               sample.
 */
template<typename P> double odara(const P &p) {
    return odara(p.cbegin(), p.cend(), p.dims());
}

} // namespace Wheelwright

#endif // WHEELWRIGHT_BOUNDS_H_
