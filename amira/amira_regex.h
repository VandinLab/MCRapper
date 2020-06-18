/**
 * Regexp used in AMIRA (support) code to parse itemsets.
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

#ifndef AMIRA_REGEX_H_
#define AMIRA_REGEX_H_

#include <regex>

namespace amira {

const static std::regex FIre {"((?:[0-9]+ )*)\\(([0-9]+)\\)",
    std::regex_constants::optimize};

const static std::regex JSONFIre {
    "\\s*\"((?:(?:[0-9]+_)*[0-9]+)|\\*)\": ([0-9]+),?",
        std::regex_constants::optimize};

const static std::regex epsilonRe { "\\s*eps: [0-9]+(.[0-9][0-9]*)?,?",
    std::regex_constants::optimize};

const static std::regex thetaRe { "\\s*minimum_frequency: [0-9]+(.[0-9][0-9]*)?,?",
    std::regex_constants::optimize};


} // namespace amira

#endif // AMIRA_REGEX_H_
