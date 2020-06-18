/**
 * Functions to compute the upper bound to the empirical Rademacher averages and
 * the maximum deviation.
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
 *
 */

#ifndef _EPSILON_H_
#define _EPSILON_H_

#include <algorithm>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <nlopt.hpp>

#include "wheelwright/matteoutils/LogSumFromLogs.h"
#include "wheelwright/matteoutils/MissingMath.h"
#include "wheelwright/Bounds.h"

#include "amira.h"

namespace amira {

// Store the bounds to the ERA and to the maximum deviation.
struct EraEps {
    double era; // deterministic upper bound to the empirical Rademacher average
    double eps; // probabilistic upper bound to the maximum deviation
};

namespace impl {

// Evaluate the objective function at x, using data as needed. The second
// argument is the gradient at x.
// This function equals the function in the paper multiplied by the sample
// size, but is for the specific case of computing the first upper
// bound to the empirical Rademacher average (also known as omega1).
inline double items_objective(const std::vector<double> &x,
        std::vector<double> &, void *data) {
    const auto &infos {*(static_cast<std::map<item,ItemsetInfo>*>(data))};
    const auto lncoshx {MatteoUtils::logcosh(x[0])};
    MatteoUtils::LogSumFromLogs sum {std::numeric_limits<double>::lowest()};
    for (const auto &info : infos) {
        const auto &h {info.second.h};
        sum += info.second.sp * lncoshx;
        for (const auto &gv : info.second.g) {
            count harminusj {h.at(gv.first)};
            for (const auto end {harminusj - gv.second}; harminusj > end;
                    --harminusj) {
                const auto twoexp {std::min(gv.first, harminusj)};
                const auto toadd {twoexp * MatteoUtils::log_2 + (1u + harminusj)
                    * lncoshx};
                sum += toadd;
            }
        }
    }
    return sum.value() / x[0];
}

// Evaluate the objective function at x, using data as needed. The second
// argument is the gradient at x.
// This function equals the function in the paper multiplied by the sample
// size.
inline double itemset_objective(const std::vector<double> &x,
        std::vector<double> &, void *data) {
    auto& [mine_supp_minus1, itmswinfos] {
        *(static_cast<std::tuple<count const &,
                std::set<amira::ItemsetWithInfo,
                amira::SuppThenInvByLengthComp> const &>*>(data))};
    const auto lncoshx {MatteoUtils::logcosh(x[0])};
    MatteoUtils::LogSumFromLogs sum {std::numeric_limits<double>::lowest()};
    // TODO: Convert to OpenMP for loop. We do the usual reduction on sum.
    for (const auto &itmswinfo : itmswinfos) {
        const auto &h {itmswinfo.info.h};
        const auto &w {itmswinfo.info.w};
        sum += itmswinfo.info.sp * lncoshx;
        for (const auto &gv : itmswinfo.info.g) {
            count harminusj {h.at(gv.first)};
            count warminusj_plus1 {w.at(gv.first) + 1};
            for (const auto end {harminusj - gv.second}; harminusj > end;
                    --harminusj, --warminusj_plus1) {
                // XXX: The following computation of the log of the sum of
                // binomials slows down everything by too much. Instead, we use
                // an upper bound from Lemma 3.8.2 of Lovász L. and Pelikán, J.
                // and Vesztergombi, K., Discrete Mathematics. The upper bound
                // only applies when harminusj is even, and 'mine_supp_minus1 -
                // 1' is at most harminus/2.
                //
                //double cis;
                //// Shortcut in case the sum of the binomial coefficents would
                //// be equal to 2^harminusj). We use ">" in the test because it
                //// is equivalent to using "mine_supp_minus1 -1" on the lhs and
                //// ">=".
                //if (mine_supp_minus1 > harminusj)
                //    cis = MatteoUtils::log_2 * std::min(gv.first, harminusj);
                //else
                //    cis = std::min(MatteoUtils::log_2 * gv.first,
                //            MatteoUtils::log_sum_of_binoms(harminusj,
                //                mine_supp_minus1 - 1, 0));
                double cis;
                if (harminusj % 2 == 0 &&
                        mine_supp_minus1 - 1 <= harminusj / 2 ) {
                    const double bound {
                        MatteoUtils::logbinom(harminusj, mine_supp_minus1 - 1)
                            - MatteoUtils::logbinom(harminusj, harminusj / 2)
                            + MatteoUtils::log_2 * (harminusj - 1)};
                    cis = std::min(MatteoUtils::log_2 * gv.first, bound);
                } else
                    cis = MatteoUtils::log_2 * std::min(gv.first,
                        harminusj);
                const auto supp {std::min(mine_supp_minus1, warminusj_plus1)};
                sum += cis + supp * lncoshx;
            }
        }
    }
    return sum.value() / x[0];
}

} // namespace impl

template<class T> void compute_eraeps(const double d, const count size,
        T &data, EraEps &res, typename std::enable_if_t<
        std::is_same<T, std::map<item,ItemsetInfo>>::value ||
        std::is_same<T, std::tuple<count const &,
                std::set<amira::ItemsetWithInfo,
                amira::SuppThenInvByLengthComp> const &>>::value,
            bool> = false) {
    // Use the infos to compute the upper bound for the Rademacher average.
    // Most of the following code comes from Wheelwright.
    // Set up the optimization problem. We use a *local* optimization algorithm
    // because the objective function, even if not convex, has a unique local
    // minimum, which is therefore also global.
    // XXX MR: Since we don't implement the gradient computation in the
    // *_objective functions, we run a derivative-free algorithm.
    nlopt::opt opt_prob(nlopt::LN_COBYLA, 1u);
    // This variable keeps the value of the smallest objective value computed
    // during the optimization. We have to keep this value because if the
    // optimization throws a nlopt::roundoff_limited exception, there is no
    // way to obtain the objective value, despite what the NLopt documentation
    // says.
    std::vector<double> lb({std::numeric_limits<double>::min()});
    opt_prob.set_lower_bounds(lb);
    opt_prob.set_ftol_abs(1e-9);
    // Set initialization point. The choice of '2' is somewhat arbitrary.
    std::vector<double> x {2.0};
    double objval {0};
    if constexpr (std::is_same<T, std::map<item,ItemsetInfo>>::value) {
        opt_prob.set_min_objective(impl::items_objective,
                &data);
    } else if constexpr (std::is_same<T,
            std::tuple<count const &,
                std::set<amira::ItemsetWithInfo,
                amira::SuppThenInvByLengthComp> const &>>::value) {
        opt_prob.set_min_objective(impl::itemset_objective,
                &data);
    }
    try {
        opt_prob.optimize(x, objval);
    } catch (const std::runtime_error &e) {
        std::stringstream ss;
        ss << "Error in the optimization: either a generic NLopt failure or "
            "the objective function computation threw a std::runtime_error. "
            "The message was: " << e.what();
        throw std::runtime_error(ss.str());
    }
    std::cout << " s opt " << x[0] << std::endl;
    // Divide by size because objective() computes the objective value
    // multiplied by the sample size.
    res.era = objval / size;
    // Computing the epsilon by first computing a (probabilistic) upper bound to
    // the Rademacher average and then use it to compute the bound to the
    // deviation.
    const auto ra {Wheelwright::ra_bound(res.era, size, d / 2, 1,
            Wheelwright::RaUBMethod::SELFBOUNDING)};
    res.eps = Wheelwright::dev_bound(ra, size, d / 2, 1,
            Wheelwright::MaxDevType::ABSOL);
}

} // namespace amira

#endif // _EPSILON_H_
