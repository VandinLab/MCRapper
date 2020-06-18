/** @file Bounds_impl.h
 * Implementation details for Bounds.h and Bounds.cpp -- nothing to see here.
 *
 * @author Matteo Riondato
 * @date 2017 07 06
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

#ifndef WHEELWRIGHT_BOUNDS_IMPL_H_
#define WHEELWRIGHT_BOUNDS_IMPL_H_

#include <cmath>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifdef WHEELWRIGHT_NLOPT
#include <nlopt.hpp>
#else
#include <tuple>
#endif

#include "matteoutils/Iterators.h"
#include "matteoutils/LogSumFromLogs.h"
#include "matteoutils/MissingMath.h"

namespace Wheelwright {

template<typename Iterator> double odara(const Iterator,
    const Iterator, const std::size_t, std::enable_if_t< !
        MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0);

template<typename Iterator> double odara(const Iterator,
    const Iterator, const std::size_t, std::enable_if_t<
        MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0);
} // namespace Wheelwright

namespace Wheelwright::impl {

/** Check whether a class has member function with signature ``double
 * elltwosquared()``.
 *
 * The implementation is inspired by the one at
 * https://stackoverflow.com/a/16824239 .
 *
 * @tparam C  the type to check.
 */
template<typename C> struct has_elltwosquared {
    // TODO: Wheelwright::impl::has_elltwosquared: We may be able to accomplish
    // this functionality differently with c++1z, using void_t. See see
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4436.pdf .
    private:
        template<typename T> static constexpr auto check(T*)
            -> typename std::is_same<double,
                decltype(std::declval<T>().elltwosquared())>::type;

        template<typename> static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:
        /** The value to use as the template argument of std::enable_if_t. */
        static constexpr bool value = type::value;
};

/** Check whether a class has member function ``domain_begin()``.
 *
 * The implementation is inspired by the one at
 * https://stackoverflow.com/a/257382 .
 *
 * @tparam C  the type to check.
 */
template<typename C> struct has_domain_begin {
    private:
        template<typename T> static constexpr char
            check(decltype(&T::domain_begin));

        template<typename> static constexpr long check(...);

    public:
        /** The value to use as the template argument of std::enable_if_t. */
        enum { value = sizeof(check<C>(0)) == sizeof(char) };
};

/** Check whether a class has member function ``ondedrawsum()``.
 *
 * The implementation is inspired by the one at
 * https://stackoverflow.com/a/257382 .
 *
 * @tparam C  the type to check.
 */
template<typename C> struct has_onedrawsum {
    private:
        template<typename T> static constexpr char
            check(decltype(&T::onedrawsum));

        template<typename> static constexpr long check(...);

    public:
        /** The value to use as the template argument of std::enable_if_t. */
        enum { value = sizeof(check<C>(0)) == sizeof(char) };
};

/** Returns the squared \f$\ell_2\f$ norm of @a v (arithmetic-type
 * specialization).
 *
 * @tparam T  the type of the object @a v. If T is an arithmetic type, returns
 *            @a v, otherwise, returns @a v.elltwosquared().
 *
 * @param v  the object whose norm to return.
 */
template<typename T> T elltwosquared(const T v,
        std::enable_if_t<std::is_arithmetic<T>::value>* = 0) {
    return v;
}

/** Returns the squared \f$\ell_2\f$-norm of @a v (non-arithmetic-type
 * specialization).
 *
 * @tparam T  the type of the object @a v. If T is an arithmetic type, returns
 *            @a v, otherwise, returns @a v.elltwosquared().
 *
 * @param v  the object whose norm to return.
 */
template<typename T> double elltwosquared(const T v,
        std::enable_if_t<has_elltwosquared<T>::value>* = 0) {
    return v.elltwosquared();
}

/** Returns the maximum \f$\ell_2\f$-norm of a sequence [@a b, @a e)
 * (non-random-access-iterator version).
 *
 * @tparam Iterator  the type of the iterators defining the sequence.
 *
 * @param b,e  Iterator objects defining the sequence.
 */
template<typename Iterator> double maxelltwosquared(Iterator b, Iterator e,
        std::enable_if_t<MatteoUtils::is_iterator<Iterator>::value>* = 0,
        std::enable_if_t<
            ! MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
    double max {0};
    for (; b != e; ++b)
        if (impl::elltwosquared(*b) > max)
            max = impl::elltwosquared(*b);
    return max;
}

/** Returns the maximum \f$\ell_2\f$-norm of a sequence [@a b, @a e)
 * (random-access-iterator version).
 *
 * @tparam Iterator  the type of the iterators defining the sequence.
 *
 * @param b,e  Iterator objects defining the sequence.
 */
template<typename Iterator> double maxelltwosquared(Iterator b, Iterator e,
        std::enable_if_t<
        MatteoUtils::is_random_access_iterator<Iterator>::value>* = 0) {
    const auto dist {std::distance(b, e)};
    double max {0};
    #pragma omp parallel for reduction(max : max)
    for (auto i = decltype(dist){0}; i < dist; ++i)
        if (impl::elltwosquared(*(b + i)) > max)
            max = impl::elltwosquared(*(b + i));
    return max;
}

template<typename T> T onedrawsum(const T v,
        std::enable_if_t<std::is_arithmetic<T>::value>* = 0) {
    return v;
}

template<typename T> double onedrawsum(const T v,
        std::enable_if_t<has_onedrawsum<T>::value>* = 0) {
    return v.onedrawsum();
}

/** Returns the value of the objective function for the AGOR and RU bounds
 * computed at @a x.
 */
inline double objective_AGORRU(const std::vector<double> &x,
        std::vector<double> &grad, void *f_data) {
#ifndef WHEELWRIGHT_NLOPT
#error "The optimization-based bounds to the empirical Rademacher average are not available because Wheelwright is being compiled without NLopt support ('WHEELWRIGHT_NLOPT' is not defined)."
    std::ignore = x;
    std::ignore = grad;
    std::ignore = f_data;
#endif // WHEELWRIGHT_NLOPT
    std::vector<double> *sum_exponents_p =
        static_cast<std::vector<double>*>(f_data);
    const double s_square {std::pow(x[0], 2)};
    MatteoUtils::LogSumFromLogs runningLogSum;
    MatteoUtils::LogSumFromLogs runningLogSumForGrad2ndTerm;
    #pragma omp declare reduction  \
        (runningLogSumReduction : MatteoUtils::LogSumFromLogs: omp_out += omp_in)
    #pragma omp parallel for \
        reduction(runningLogSumReduction : runningLogSum, \
                runningLogSumForGrad2ndTerm)
    for (std::size_t i = 0; i < (*sum_exponents_p).size(); ++i) {
        const double s_expo {s_square * (*sum_exponents_p)[i]};
        runningLogSum += s_expo;
        if (! grad.empty())
            runningLogSumForGrad2ndTerm += s_expo + std::log((*sum_exponents_p)[i]);
    }
    const double logSum {runningLogSum.value()};
    if (! grad.empty())
        grad[0] = 2 * std::exp(runningLogSumForGrad2ndTerm.value() - logSum)
            - logSum / s_square;
    return logSum / x[0];
}

inline double objective_AGORRU_FD(const std::vector<double> &x,
        std::vector<double> &grad, void *f_data) {
    // If x is very very small, just return a big number and a fast decreasing
    // gradient.
    if (x[0] < 1e-300) {
        if (! grad.empty())
            grad[0] = -100;
        return std::numeric_limits<double>::max();
    }

    const auto *data {static_cast<std::pair<std::vector<double>,
                std::vector<std::vector<std::size_t>>>*>(f_data)};
    const auto &args {data->first};
    const auto &counts {data->second};
    // XXX It would be nice to use a binding as below, but because we use counts
    // in the OpenMP loop later, counts would be captured in a lambda, and that
    // is not allowed for bindings. See also
    // https://stackoverflow.com/a/46115028/.
    //const auto & [args, counts] =
    //    *(static_cast<std::pair<std::vector<double>,
    //            std::vector<std::vector<std::size_t>>>*>(f_data));

    std::vector<double> log_cosh_x_by_args;
    std::vector<double> log_args_by_sinh_x_by_args;
    for (const auto & arg : args) {
        const auto x_by_arg {x.at(0) * arg};
        log_cosh_x_by_args.emplace_back(MatteoUtils::logcosh(x_by_arg));
        log_args_by_sinh_x_by_args.emplace_back(
                std::log(arg * std::sinh(x_by_arg)));
    }

    double funcLS_term {0};
    auto count_it {counts.at(0).begin()};
    auto log_cosh_x_by_args_it {log_cosh_x_by_args.begin()};
    while (count_it != counts.at(0).end())
        funcLS_term += *(log_cosh_x_by_args_it++) * *(count_it++);
    MatteoUtils::LogSumFromLogs funcLS {funcLS_term};

    MatteoUtils::LogSumFromLogs gradLS;
    std::size_t start_idx {0};
    if (! grad.empty()) {
        if (funcLS_term == 0) {
            // We cannot initialize the gradient with the first vector because
            // it is all zeroes, so find the first one that isn't.
            // XXX Although really it should just be the next one because we are
            // assuming to be working with a *set* of vectors.
            // XXX We are also assuming not all vectors are zero.
            for (start_idx = 1; true; ++start_idx) {
                funcLS_term = 0;
                count_it = counts.at(start_idx).begin();
                log_cosh_x_by_args_it = log_cosh_x_by_args.begin();
                while (count_it != counts.at(start_idx).end())
                    funcLS_term += *(log_cosh_x_by_args_it++) * *(count_it++);
                funcLS += funcLS_term;
                if (funcLS_term != 0)
                    break;
            }
        }
        auto log_args_by_sinh_x_by_args_it {log_args_by_sinh_x_by_args.begin()};
        count_it = counts.at(start_idx).begin();
        log_cosh_x_by_args_it = log_cosh_x_by_args.begin();

        // Find the first positive count, which must exists.
        int to_skip = 0;
        while (*count_it == 0) {
            ++count_it;
            ++log_args_by_sinh_x_by_args_it;
            ++log_cosh_x_by_args_it;
            ++to_skip;
        }
        auto count_it_prev {count_it};
        auto log_cosh_x_by_args_it_prev {log_cosh_x_by_args_it};
        auto log_args_by_sinh_x_by_args_it_prev {log_args_by_sinh_x_by_args_it};
        auto gradLS_term_term {funcLS_term - *(log_cosh_x_by_args_it++) +
            *(log_args_by_sinh_x_by_args_it++) + std::log(*(count_it++))};
        MatteoUtils::LogSumFromLogs gradLS_term {gradLS_term_term};
        while (count_it != counts.at(start_idx).end()) {
            if (*count_it != 0) {
                gradLS_term_term = gradLS_term_term
                    - *log_args_by_sinh_x_by_args_it_prev
                    - std::log(*count_it_prev)
                    + *log_cosh_x_by_args_it_prev - *log_cosh_x_by_args_it
                    + *log_args_by_sinh_x_by_args_it
                    + std::log(*count_it);
                gradLS_term += gradLS_term_term;
                // Only increment the *_prev iterators when count_it != 0,
                // otherwise at the next iteration we take the logarithm of
                // zero, which is not good.
                count_it_prev = count_it++;
                log_cosh_x_by_args_it_prev = log_cosh_x_by_args_it++;
                log_args_by_sinh_x_by_args_it_prev =
                    log_args_by_sinh_x_by_args_it++;
            } else {
                ++count_it;
                ++log_cosh_x_by_args_it;
                ++log_args_by_sinh_x_by_args_it;
            }
        }
        gradLS.set(gradLS_term.value());
    } // end if (!grad.emtpy())

    #pragma omp declare reduction (runningLogSumReduction: MatteoUtils::LogSumFromLogs: omp_out += omp_in)
    #pragma omp parallel for reduction(runningLogSumReduction:funcLS, gradLS)
    for (std::size_t i = start_idx + 1; i < counts.size(); ++i) {
        double funcLS_term_r {0};
        auto count_it_r {counts.at(i).begin()};
        auto log_cosh_x_by_args_it_r {log_cosh_x_by_args.begin()};
        while (count_it_r != counts.at(i).end())
            funcLS_term_r += *(log_cosh_x_by_args_it_r++) * *(count_it_r++);
        funcLS += funcLS_term_r;
        // We don't update the grad if it is a vector of all zeroes.
        if (funcLS_term_r > 0 && ! grad.empty()) {
            auto count_it_r {counts.at(i).begin()};
            auto log_cosh_x_by_args_it_r {log_cosh_x_by_args.begin()};
            auto log_args_by_sinh_x_by_args_it_r {
                log_args_by_sinh_x_by_args.begin()};

            int to_skip_r {0};
            // Find the first positive counter. There must be one or we would
            // have skipped this vector.
            while (*count_it_r == 0) {
                ++count_it_r;
                ++log_cosh_x_by_args_it_r;
                ++log_args_by_sinh_x_by_args_it_r;
                ++to_skip_r;
            }
            auto count_it_prev_r {count_it_r};
            auto log_cosh_x_by_args_it_prev_r {log_cosh_x_by_args_it_r};
            auto log_args_by_sinh_x_by_args_it_prev_r {
                log_args_by_sinh_x_by_args_it_r};

            auto gradLS_term_term_r {funcLS_term_r -
                *(log_cosh_x_by_args_it_r++) +
                    *(log_args_by_sinh_x_by_args_it_r++) +
                    std::log(*(count_it_r++))};
            MatteoUtils::LogSumFromLogs gradLS_term_r {gradLS_term_term_r};

            while(count_it_r != counts.at(i).end()) {
                if (*count_it_r != 0) {
                    gradLS_term_term_r = gradLS_term_term_r
                        - *log_args_by_sinh_x_by_args_it_prev_r
                        - std::log(*count_it_prev_r)
                        + *log_cosh_x_by_args_it_prev_r
                        - *log_cosh_x_by_args_it_r
                        + *log_args_by_sinh_x_by_args_it_r
                        + std::log(*count_it_r);
                    gradLS_term_r += gradLS_term_term_r;
                    // Only increment the *_prev_r iterators when count_it != 0,
                    // otherwise at the next iteration we take the logarithm of
                    // zero, which is not good.
                    count_it_prev_r = count_it_r++;
                    log_cosh_x_by_args_it_prev_r = log_cosh_x_by_args_it_r++;
                    log_args_by_sinh_x_by_args_it_prev_r =
                        log_args_by_sinh_x_by_args_it_r++;
                } else {
                    ++count_it_r;
                    ++log_cosh_x_by_args_it_r;
                    ++log_args_by_sinh_x_by_args_it_r;
                }
            }
            gradLS += gradLS_term_r.value();
        }
    }
    const double val {funcLS.value() / x[0]};
    if (! grad.empty()) {
        grad[0] = std::exp(gradLS.value() - std::log(x[0]) - funcLS.value()) -
            val / x[0];
    }
    return val;
}

/** Setup and solve the optimization problem to compute upper bounds to the
 * empirical Rademacher average.
 *
 * This function is available only when the `WHEELWRIGHT_NLOPT` preprocessor
 * directive is defined.
 *
 * @param objective  a pointer to the objective function.
 * @param data       a pointer to the data to be passed to the objective
 *                   function.
 * @param ftol       the absolute tolerance in the objective value for early
 *                   termination of the optimization.
 */
inline double optimize(double (&objective)(
            const std::vector<double> &, std::vector<double> &, void *),
        void *data, const double ftol = 1e-7) {
#ifndef WHEELWRIGHT_NLOPT
#error "The optimization-based bounds to the empirical Rademacher average are not available because Wheelwright is being compiled without NLopt support ('WHEELWRIGHT_NLOPT' is not defined)."
    std::ignore = objective;
    std::ignore = data;
    std::ignore = ftol;
#endif // WHEELWRIGHT_NLOPT
    // Use a *local* optimization algorithm because the functions, even if not
    // convex, have a unique local minimum, which is therefore also global.
    nlopt::opt opt_prob(nlopt::LD_MMA, 1u);
    opt_prob.set_min_objective(objective, data);
    std::vector<double> lb({std::numeric_limits<double>::min()});
    opt_prob.set_lower_bounds(lb);
    // The upper bound is needed for some algos, e.g. GN_DIRECT_L
    //std::vector<double> ub {std::numeric_limits<double>::max()};
    //opt_prob.set_upper_bounds(ub);
    //opt_prob.set_ftol_rel(0.01); // one percent relative tolerance
    if (ftol > 0)
        opt_prob.set_ftol_abs(ftol);
    // Set initialization point. The choice of '2' is somewhat arbitrary: the
    // optimization converges fast in our tests.
    std::vector<double> x {2.0};
    double res {0};
    try {
        opt_prob.optimize(x, res);
    } catch (std::runtime_error &e) {
        std::stringstream ss;
        ss << "Wheelwright::impl::optimize: Error in the optimization: either "
            "a generic NLopt failure or the objective function computation "
            "threw a std::runtime_error. The message was: " << e.what();
        throw std::runtime_error(ss.str());
    }
#ifdef WHEELWRIGHT_OPT_DEBUG
    // The following code tests that the optimal solution returned by NLopt is
    // indeed optimal, and that the computation of the gradient is correct.
    //
    // This code could be somehow moved to to a test.
    double step {x[0] / 10};
    for (std::size_t i {1}; i < 10; ++i) {
        std::vector<double> y {i * step};
        std::vector<double> gy {0};
        auto r {objective(y, gy, data)};
        if (r < res) {
            std::stringstream ss;
            ss << "Wheelwright::impl::optimize: the NLopt \"optimal\" "
                "objective value " << res << " for x=" << x[0] << " is not "
                "really optimal: found objective value " << r << " for x="
                << y[0] << ".";
            throw std::runtime_error(ss.str());
        }
        if (gy[0] > 0) {
            std::stringstream ss;
            ss << "Wheelwright::impl::optimize: the NLopt gradient is not "
                "decreasing as expected: found positive gradient " << gy[0]
                << " at " << y[0] << " (\"optimal\" objective value: " << res
                << " at " << x[0] << ").";
            throw std::runtime_error(ss.str());
        }
        std::vector<double> yDiff(1);
        std::vector<double> empty;
        const double h {std::sqrt(std::numeric_limits<double>::epsilon()) * y[0]};
        volatile double yph {y[0] + h};
        const double dy {yph - y[0]};
        yDiff[0] = yph;
        const double rPlus {objective(yDiff, empty, data)};
        volatile double ymh {y[0] - h};
        yDiff[0] = ymh;
        const double rMinus {objective(yDiff, empty, data)};
        const double finiteDiffQuot {(rPlus - rMinus) / (2 * dy)};
        // Arbitrarily say that the gradient is computed correctly if the
        // finiteDiffQuot is within 5% of it and the absolute difference is not
        // greater than 0.1. These are empirical errors, so sometimes the check
        // may fail.
        // The comparison between gy[0] and the finiteDiffQuot is meaningful
        // only if the objective value is greater than one of the two
        // "perturbed" objective values.
        if (r > std::min(rPlus, rMinus) && \
                std::fabs(gy[0] - finiteDiffQuot) / std::fabs(gy[0]) \
                > 0.05 && std::fabs(gy[0] - finiteDiffQuot) > 0.1)
            throw std::runtime_error("Wheelwright::impl::optimize: the "
                    "gradient computation does not appear to be correct.");
    }
    step = x[0];
    for (std::size_t i {2}; i <= 11; ++i) {
        std::vector<double> y {i * step};
        std::vector<double> gy {0};
        auto r {objective(y, gy, data)};
        if (r < res) {
            std::stringstream ss;
            ss << "Wheelwright::impl::optimize: the NLopt \"optimal\" "
                "objective value " << res << " for x=" << x[0] << " is not "
                "really optimal: found objective value " << r << " for x="
                << y[0] << ".";
            throw std::runtime_error(ss.str());
        }
        if (gy[0] < 0) {
            std::stringstream ss;
            ss << "Wheelwright::impl::optimize: the NLopt gradient is not "
                "decreasing as expected: found negative gradient " << gy[0]
                << " at " << y[0] << " (\"optimal\" objective value: " << res
                << " at " << x[0] << ").";
            throw std::runtime_error(ss.str());
        }
        std::vector<double> yDiff(1);
        std::vector<double> empty;
        const double h {std::sqrt(std::numeric_limits<double>::epsilon()) * y[0]};
        volatile double yph {y[0] + h};
        const double dy {yph - y[0]};
        yDiff[0] = yph;
        const double rPlus {objective(yDiff, empty, data)};
        volatile double ymh {y[0] - h};
        yDiff[0] = ymh;
        const double rMinus {objective(yDiff, empty, data)};
        const double finiteDiffQuot {(rPlus - rMinus) / (2 * dy)};
        // Arbitrarily say that the gradient is computed correctly if the
        // finiteDiffQuot is within 5% of it and the absolute difference is not
        // greater than 0.1. These are empirical errors, so sometimes the check
        // may fail.
        // The comparison between gy[0] and the finiteDiffQuot is meaningful
        // only if the objective value is greater than one of the two
        // "perturbed" objective values.
        if (r > std::min(rPlus, rMinus) && \
                std::fabs(gy[0] - finiteDiffQuot) / std::fabs(gy[0]) \
                > 0.05 && std::fabs(gy[0] - finiteDiffQuot) > 0.1)
            throw std::runtime_error("Wheelwright::impl::optimize: the "
                    "gradient computation does not appear to be correct.");
    }
#endif // WHEELWRIGHT_OPT_DEBUG
    return res;
}

/** Returns an upper bound to the Rademacher average using the empirical
 * Rademacher average or the odara and McDiarmid's bounded difference
 * inequality.
 *
 * @param c     the empirical Rademacher average (or an upper bound) or the
 *              odara.
 * @param size  the sample size.
 * @param eta   the probability of error (i.e., 1 - confidence).
 * @param span  the width of the co-domain of the functions (or an upper bound).
 */
inline double ra_bound_MCDIARMIDONEDRAW(const double c, const std::size_t size,
        const double eta, const double span) {
    return c + span * std::sqrt(- std::log(eta) / (2 * size));
}

/** Returns an upper bound to the Rademacher average using the empirical
 * Rademacher average and the tail bound for self-bounding functions.
 *
 * @param era   the empirical Rademacher average (or an upper bound).
 * @param size  the sample size.
 * @param eta   the probability of error (i.e., 1 - confidence).
 * @param span  the width of the co-domain of the functions (or an upper bound).
 */
inline double ra_bound_SELFBOUNDING(const double era, const std::size_t size,
        const double eta, const double span) {
    const double gamma {std::log(1 / eta)};
    const double span_gamma {span * gamma};
    const double second_term {(span_gamma + std::sqrt(span_gamma *
                (span_gamma + 4 * size * era))) / (2 * size)};
    return era + second_term;
}

} // namespace Wheelwright::impl

#endif // WHEELWRIGHT_BOUNDS_IMPL_H_
