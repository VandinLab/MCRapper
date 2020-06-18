/** @file Sampler.h
 * Generate uniformly distributed samples from 0 to n-1.
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

#ifndef MATTEOUTILS_SAMPLER_H_
#define MATTEOUTILS_SAMPLER_H_

#include <cstddef>
#include <random>
#include <vector>

namespace MatteoUtils {

/** Generate uniformly distributed independent samples with replacement from
 * 0 to @a n -1.
 *
 * The upper limit of the domain can be changed during the lifetime of the
 * object.
 */
class Sampler {
    public:
        /** Constructor.
         *
         * @param n  the size of the population. The population will be the set of
         *         integers from 0 to @a n-1.
         */
        explicit Sampler(const std::size_t n) : prng(std::random_device {}()),
                 dist(0, n - 1) {
            // TODO Sampler::Sampler : Look at
            // http://www.pcg-random.org/posts/cpp-seeding-surprises.html
            // to fix the seeding of the PRNG.
#ifdef MATTEOUTILS_INVARIANTS
            if (n < 1)
                throw std::out_of_range("MatteoUtils::Sampler::Sampler: size "
                        "must be positive");
#endif
        }

        /** Get @a n samples.
         *
         * @param n  the desired number of samples.
         *
         * @returns  a vector of @a n non-negative unsigned integers, sampled
         *           independently and uniformly at random with replacement from
         *           the population.
         */
        std::vector<std::size_t> samples(const std::size_t n) {
            std::vector<std::size_t> indices(n);
            // XXX: Sampler::samples: we cannot use OpenMP here
            // because the PRNG is not thread-safe.
            //#pragma omp parallel for
            for (std::size_t i = 0; i < n; ++i)
                indices[i] = dist(prng);
            return indices;
        }

        /** Set the random seed.
         *
         * @param s the random seed.
         */
        void seed(const unsigned int s) {
            // TODO: Sampler::seed: Look at
            // http://www.pcg-random.org/posts/cpp-seeding-surprises.html
            // to fix the seeding of the PRNG.
            prng.seed(s);
        }

        /** Set the size of the population
         *
         * @param n  the new size of the population, which will be the set of
         *         integers from 0 to @a n-1.
         */
        void set_size(const std::size_t n) {
#ifdef MATTEOUTILS_INVARIANTS
            if (n < 1)
                throw std::out_of_range("MatteoUtils::Sampler::set_size: size "
                        "must be positive");
#endif
            std::uniform_int_distribution<std::size_t> dist2(0, n - 1);
            dist.param(dist2.param());
        }

    private:
        std::mt19937_64 prng;
        std::uniform_int_distribution<std::size_t> dist;

}; // class Sample

}// namespace MatteoUtils

#endif // MATTEOUTILS_SAMPLER_H_
