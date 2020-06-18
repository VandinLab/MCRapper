/**
* Compare a collection A of itemsets with a "ground-truth" collection B of
* itemsets, to check that A is an eps-approximation of B, for eps passed in the
* input of the collection A. The input to this utility must be the output of
* amira or grahne/fim_*, passed through the sort_fis utility.
*
* Copyright 2015, 2018 Matteo Riondato <riondato@acm.org>
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

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <unistd.h>

#include "amira.h"
#include "amira_regex.h"

void usage(const char *binary) {
    std::cerr << binary << ": compare two collections of frequent itemsets"
        << std::endl
        << "USAGE: " << binary << " [-hjv] groundtruth approx" << std::endl
        << "\t-h : print this message and exit" << std::endl
        << "\t-j : approx input is in JSON format" << std::endl
        << "\t-J : groundtruth input is in JSON format" << std::endl
        << "\t-o : output in JSON format" << std::endl
        << "\t-v : print log messages to stderr during the execution"
        << std::endl;
}

int main(int argc, char **argv) {
    const auto start {std::chrono::system_clock::now()};
    bool verbose {false};
    bool json_app {false};
    bool json_gro {false};
    bool json_out {false};
    {
        extern int optind;
        char opt;
        while ((opt = getopt(argc, argv, "hjJov")) != -1) {
            switch (opt) {
                case 'h':
                    usage(argv[0]);
                    return EXIT_SUCCESS;
                    break;
                case 'j':
                    json_app = true;
                    break;
                case 'J':
                    json_gro = true;
                    break;
                case 'o':
                    json_out = true;
                    break;
                case 'v':
                    verbose = true;
                    break;
                default:
                    std::cerr << "Error: wrong option." << std::endl;
                    return EXIT_FAILURE;
            }
        }
        if (optind != argc - 2) {
            std::cerr << "Error: wrong number of arguments" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::ifstream g {argv[argc - 2]};
    if (g.fail()) {
        std::cerr << "Error: can't open the groundtruth file '"
            << argv[argc - 2] << "'." << std::endl;
        return EXIT_FAILURE;
    }
    std::ifstream a {argv[argc - 1]};
    if (a.fail()) {
        std::cerr << "Error: can't open the approximation file '"
            << argv[argc - 1] << "'." << std::endl;
        return EXIT_FAILURE;
    }
    // Read the approximation, finding epsilon and theta, and reading all the
    // itemsets with frequency at least theta - eps.
    if (verbose)
        std::cerr << "Looking for theta...";
    std::string line;
    std::cmatch m;
    // Skip the heading invalid lines until we find the theta line.
    while (std::getline(a, line) &&
            ! std::regex_match(line.c_str(), m, amira::thetaRe))
        ; // empty body loop
    if (a.fail()) {
        std::cerr << "Error: malformed input: can't find theta "
            "(i.e., minimum_frequency) in the approximation file." << std::endl;
        return EXIT_FAILURE;
    }
    // line now contain the theta line. Process it and get the theta.
    const double theta {std::strtod(m[1].first, NULL)};
    if (errno == ERANGE || theta <= 0) {
        std::cerr << "Error: malformed input: theta must be positive."
            << std::endl;
        return EXIT_FAILURE;
    }
    if (verbose)
        std::cerr << "done: " << theta << std::endl;
    // Skip the heading invalid lines until we find the epsilon line.
    if (verbose)
        std::cerr << "Looking for epsilon...";
    while (std::getline(a, line) &&
            ! std::regex_match(line.c_str(), m, amira::epsilonRe))
        ; // empty body loop
    if (a.fail()) {
        std::cerr << "Error: malformed input: can't find epsilon in the "
            "approximation file." << std::endl;
        return EXIT_FAILURE;
    }
    // line now contain the epsilon line. Process it and get the epsilon.
    const double eps {std::strtod(m[1].first, NULL)};
    if (errno == ERANGE || eps <= 0) {
        std::cerr << "Error: malformed input: epsilon must be positive."
            << std::endl;
        return EXIT_FAILURE;
    }
    if (verbose)
        std::cerr << "done: " << eps << std::endl;
    if (json_app) {
        if (verbose)
            std::cerr << "Looking for itemsets in JSON approximation input...";
        // Skip the lines until we find a valid itemset line
        while (std::getline(a, line) &&
                ! std::regex_match(line.c_str(), m, amira::JSONFIre))
            ; // empty body loop
    } else { // not JSON input
        if (verbose)
            std::cerr << "Looking for itemsets in plaintext approximation "
                "input...";
        // Skip the lines until we find a valid itemset line
        while (std::getline(a, line) &&
                ! std::regex_match(line.c_str(), m, amira::FIre))
            ; // empty body loop
    }
    if (a.fail()) {
        if (verbose)
            std::cerr << std::endl;
        std::cerr << "Error: malformed input: can't find any itemset in the "
            "approximation file." << std::endl;
        return EXIT_FAILURE;
    }
    // The size of the sample.
    const amira::count asz {std::strtoul(m[2].first, NULL, 10)};
    if (errno == ERANGE || asz == 0) {
        if (verbose)
            std::cerr << std::endl;
        std::cerr << "Error: malformed input: the size of the sample must be "
            "positive." << std::endl;
        return EXIT_FAILURE;
    }
    double lb {theta - eps};
    if (lb < 0)
        lb = 0;
    std::unordered_map<std::string,double> af;
    if (json_app) {
        while (std::getline(a, line) &&
                std::regex_match(line.c_str(), m, amira::JSONFIre)) {
            const auto s {std::strtoul(m[2].first, NULL, 10)};
            const auto f {static_cast<double>(s) / asz};
            if (f >= lb)
                af.emplace(m[1], f);
            else
                break;
        }
    } else { // not JSON input
        while (std::getline(a, line) &&
                std::regex_match(line.c_str(), m, amira::FIre)) {
            const auto s {std::strtoul(m[2].first, NULL, 10)};
            const auto f {static_cast<double>(s) / asz};
            if (f >= lb) {
                std::string i {m[1]};
                i.pop_back(); // remove the last space;
                std::replace(i.begin(), i.end(), ' ', '_');
                af.emplace(i, f);
            } else
                break;
        }
    }
    if (af.empty()) {
        if (verbose)
            std::cerr << std::endl;
        std::cerr << "Error: malformed input: can't find any itemset in the "
            "approximation file." << std::endl;
        return EXIT_FAILURE;
    }
    if (verbose)
        std::cerr << "done: " << af.size() << " itemsets in approximation"
            << std::endl;
    a.close();
    if (json_gro) {
        if (verbose)
            std::cerr << "Comparing with itemsets in JSON groundtruth input...";
        // Skip the lines until we find a valid itemset line, which is the
        // emptyset (in the normal case of comparing FIs, not CFIs).
        while (std::getline(g, line) &&
                ! std::regex_match(line.c_str(), m, amira::JSONFIre))
            ; // empty body loop
    } else { // not JSON input
        if (verbose)
            std::cerr << "Comparing with itemsets in plaintext groundtruth "
                "input...";
        // Skip the lines until we find a valid itemset line, which is the
        // emptyset (in the normal case of comparing FIs, not CFIs).
        while (std::getline(g, line) &&
                ! std::regex_match(line.c_str(), m, amira::FIre))
            ; // empty body loop
    }
    if (g.fail()) {
        if (verbose)
            std::cerr << std::endl;
        std::cerr << "Error: malformed input: can't find any itemset in "
            "the groundtruth file." << std::endl;
        return EXIT_FAILURE;
    }
    // At this point the match contains the size of the original dataset.
    const amira::count gsz {std::strtoul(m[2].first, NULL, 10)};
    // Read the ground truth. We do not need to store the itemsets in it.
    // Rather, because the list is sorted, we can check them one by one wrt to
    // those in the approximation: each time we find one, we can remove it from
    // the approximation. Before that, we should check that the error in the
    // frequency is at most epsilon (or whatever it should be). Once the one we
    // read from the ground truth has frequency less than theta, the remainder
    // are false positives. We still need to continue the check and remove until
    // we read from the ground truth an itemsets with frequency less than
    // theta-epsilon.
    unsigned int orig {0}; // frequent itemsets in the ground truth
    unsigned int tp {0}; // true positives
    unsigned int fp {0}; // false positives
    unsigned int possfp {0}; // possible false positives
    decltype(af) fn; // false negatives
    std::vector<double> errs; // frequency errors for the true positives and
                              // acceptables false positives
    double err_sum {0}; // sum of the errors.
    unsigned int wrongeps {0}; // errors greater than eps
    if (json_gro) {
        // Check the true positives (i.e., real FIs in the ground truth.)
        while (std::getline(g, line) &&
                std::regex_match(line.c_str(), m, amira::JSONFIre)) {
            const amira::count s {std::strtoul(m[2].first, NULL, 10)};
            const double f {static_cast<double>(s) / gsz};
            if (f < theta)
                break;
            ++orig;
            if (const auto afn {af.extract(m[1])}; ! afn.empty()) {
                ++tp;
                const auto err {std::abs(afn.mapped() - f)};
                if (err > eps / 2)
                    ++wrongeps;
                err_sum += err;
                errs.emplace_back(err);
            } else
                fn.emplace(m[1], s);
        }
        // Check the acceptable false positives (i.e., itemsets in the ground
        // truth with frequency at most lb = theta - eps.)
        if (std::regex_match(line.c_str(), m, amira::JSONFIre)) {
            // We use a do...while because the last line read was valid, but not
            // processed.
            do {
                ++possfp;
                const amira::count s {std::strtoul(m[2].first, NULL, 10)};
                const double f {static_cast<double>(s) / gsz};
                if (f < lb)
                    break;
                if (const auto afn {af.extract(m[1])}; ! afn.empty()) {
                    ++fp;
                    const auto err {std::abs(afn.mapped() - f)};
                    if (err > eps / 2)
                        ++wrongeps;
                    err_sum += err;
                    errs.emplace_back(err);
                }
            } while (std::getline(g, line) &&
                    std::regex_match(line.c_str(), m, amira::JSONFIre));
        }
        if (! std::regex_match(line.c_str(), m, amira::JSONFIre)) {
            if (verbose)
                std::cerr << std::endl;
            std::cerr << "Warning: reached the end of the itemsets in the "
                "ground truth possibly *before* having finished to read "
                "all the possible acceptable false positives." << std::endl;
        }
    } else { // not JSON input
        // Check the true positives (i.e., real FIs in the ground truth.)
        while (std::getline(g, line) &&
                std::regex_match(line.c_str(), m, amira::FIre)) {
            const amira::count s {std::strtoul(m[2].first, NULL, 10)};
            const double f {static_cast<double>(s) / gsz};
            if (f < theta)
                break;
            ++orig;
            std::string i {m[1]};
            i.pop_back(); // remove last space
            std::replace(i.begin(), i.end(), ' ', '_');
            if (auto afn {af.extract(i)}; ! afn.empty()) {
                ++tp;
                const auto err {std::abs(afn.mapped() - f)};
                if (err > eps / 2)
                    ++wrongeps;
                err_sum += err;
                errs.emplace_back(err);
            } else
                fn.emplace(i, s);
        }
        // Check the acceptable false positives (i.e., itemsets in the ground
        // truth with frequency at most lb = theta - eps.)
        if (std::regex_match(line.c_str(), m, amira::FIre)) {
            // We use a do...while because the last line read was valid, but not
            // processed.
            do {
                ++possfp;
                const amira::count s {std::strtoul(m[2].first, NULL, 10)};
                const double f {static_cast<double>(s) / gsz};
                if (f < lb)
                    break;
                std::string i {m[1]};
                i.pop_back(); // remove last space
                std::replace(i.begin(), i.end(), ' ', '_');
                if (const auto afn {af.extract(i)}; ! afn.empty()) {
                    ++fp;
                    const auto err {std::abs(afn.mapped() - f)};
                    if (err > eps)
                        ++wrongeps;
                    err_sum += err;
                    errs.emplace_back(err);
                }
            } while (std::getline(g, line) &&
                    std::regex_match(line.c_str(), m, amira::FIre));
        }
        if (! std::regex_match(line.c_str(), m, amira::FIre)) {
            if (verbose)
                std::cerr << std::endl;
            std::cerr << "Warning: reached the end of the itemsets in the "
                "ground truth possibly *before* having finished to read "
                "all the possible acceptable false positives." << std::endl;
        }
    }
    if (verbose)
        std::cerr << "done (tp: " << tp << ", fn: " << fn.size() << ", fp: "
            << fp;
    // Anything that is left in af at this point is a non-acceptable false
    // positive.
    const auto nafp {af.size()};
    // Fraction of total possible acceptable false positives in the approx;
    const auto fpfrac {static_cast<double>(fp) / possfp};
    // Error distribution statistics
    const auto err_avg {static_cast<double>(err_sum) / errs.size()};
    const auto err_centrsqsum {std::accumulate(errs.cbegin(), errs.cend(),
            // Using '0.0' is necessary to have the return type being double.
            0.0, [err_avg] (const double a, const double b) -> auto {
            return a + std::pow(b - err_avg, 2);})};
    // Population standard deviation
    const auto err_stdev {std::sqrt(err_centrsqsum / errs.size())};
    // Sorting allows us to compute the quartiles.
    std::sort(errs.begin(), errs.end());
    const auto err_min {errs.front()};
    const auto err_max {errs.back()};
    double err_med {errs[errs.size() / 2]};
    if (errs.size() % 2 == 0)
        err_med = (err_med + errs[errs.size() / 2 - 1]) / 2;
    // XXX MR: This quartile computation is a bit rough, and may not be entirely
    // correct, but it should be good enough.
    const auto frstqidx {(errs.size() / 2 + 1) / 2};
    double err_1q {errs[frstqidx]};
    if ((errs.size() / 2 + 1) % 2 == 0)
        err_1q = (err_1q + errs[frstqidx - 1]) / 2;
    const auto thrdqidx {(errs.size() / 2) + frstqidx};
    double err_3q {errs[thrdqidx]};
    if ((errs.size() / 2 + frstqidx) % 2 == 0)
        err_3q = (err_3q + errs[frstqidx + 1]) / 2;
    if (verbose)
        std::cerr << ", nafp: " << nafp << ", wrongeps: " << wrongeps << ")"
            << std::endl << "Printing output and exiting. Goodbye."
            << std::endl;
    // Print results.
    std::string_view comma;
    std::string_view quotes;
    std::string_view sep {" "};
    std::string_view open_supp {"("};
    std::string_view close_supp {")"};
    std::string_view tab;
    if (json_out) {
        comma = ",";
        quotes = "\"";
        sep = "_";
        open_supp = ": ";
        close_supp = "";
        tab = "\t";
        std::cout << "{" << std::endl; // open everything
    } else
        std::cout << "## Compare FIs" << std::endl << std::endl;
    const auto starttime {std::chrono::system_clock::to_time_t(start)};
    std::cout << tab << quotes << "date" << quotes << ": " << quotes
        << std::put_time(std::localtime(&starttime), "%F %T") << quotes
        << comma << std::endl;
    std::cout << tab << quotes << "groundtruth" << quotes << ": " << quotes
        << argv[argc -2] << quotes << comma << std::endl;
    std::cout << tab << quotes << "approximation" << quotes << ": " << quotes
        << argv[argc -1] << quotes << comma << std::endl;
    std::cout << tab << quotes << "theta" << quotes << ": " << theta << comma
        << std::endl;
    std::cout << tab << quotes << "epsilon" << quotes << ": " << eps << comma
        << std::endl;
    std::cout << tab << quotes << "lb" << quotes << ": " << lb << comma
        << std::endl;
    std::cout << tab << quotes << "orig" << quotes << ": " << orig << comma
        << std::endl;
    std::cout << tab << quotes << "tp" << quotes << ": " << tp << comma
        << std::endl;
    std::cout << tab << quotes << "fn" << quotes << ": " << fn.size() << comma
        << std::endl;
    std::cout << tab << quotes << "fp" << quotes << ": " << fp << comma
        << std::endl;
    std::cout << tab << quotes << "possfp" << quotes << ": " << possfp << comma
        << std::endl;
    std::cout << tab << quotes << "nafp" << quotes << ": " << nafp << comma
        << std::endl;
    std::cout << tab << quotes << "fpfrac" << quotes << ": " << fpfrac << comma
        << std::endl;
    std::cout << tab << quotes << "wrongeps" << quotes << ": " << wrongeps
        << comma << std::endl;
    std::cout << tab << quotes << "err_avg" << quotes << ": " << err_avg
        << comma << std::endl;
    std::cout << tab << quotes << "err_stdev" << quotes << ": " << err_stdev
        << comma << std::endl;
    std::cout << tab << quotes << "err_min" << quotes << ": " << err_min
        << comma << std::endl;
    std::cout << tab << quotes << "err_max" << quotes << ": " << err_max
        << comma << std::endl;
    std::cout << tab << quotes << "err_med" << quotes << ": " << err_med
        << comma << std::endl;
    std::cout << tab << quotes << "err_1q" << quotes << ": " << err_1q << comma
        << std::endl;
    std::cout << tab << quotes << "err_3q" << quotes << ": " << err_3q;
    // Print the false negatives and the non-acceptable false positives.
    if (! fn.empty()) {
        if (json_out)
            std::cout << comma << std::endl << tab << quotes << "fns"
                << quotes << ": {" << std::endl;
        else
            std::cout << std::endl << std::endl << "### False Negatives"
                << std::endl << std::endl;
        // Iterate over the elements of fn to print them, but stop before the
        // last one
        auto it {fn.cbegin()};
        for (std::size_t proc {1}; proc < fn.size(); ++proc, ++it) {
            const double f {static_cast<double>(it->second) / gsz};
            std::cout << tab << tab << quotes << it->first << quotes
                << open_supp << f << close_supp << comma << std::endl;
        }
        // Print the last element of fn
        const double f {static_cast<double>(it->second) / gsz};
        std::cout << tab << tab << quotes << it->first << quotes
            << open_supp << f << close_supp << std::endl;
        if (json_out) {
            std::cout << tab <<  "}";
            if (af.empty())
                std::cout << std::endl;
        }
    }
    if (! af.empty()) {
        std::cout << comma << std::endl;
        if (json_out)
            std::cout << tab << quotes << "nafps" << quotes << ": {"
                << std::endl;
        else {
            if (fn.empty())
                std::cout << std::endl;
            std::cout << "### Non-acceptable False Positives" << std::endl
                << std::endl;
        }
        // Iterate over the elements of af to print them, but stop before the
        // last one
        auto it {af.cbegin()};
        for (std::size_t proc {1}; proc < af.size(); ++proc, ++it)
            std::cout << tab << tab << quotes << it->first << quotes
                << open_supp << it->second << close_supp << comma << std::endl;
        // Print the last element of af
        std::cout << tab << tab << quotes << it->first << quotes
            << open_supp << it->second << close_supp << std::endl;
        if (json_out)
            std::cout << tab <<  "}" << std::endl;
    }
    if (fn.empty() && af.empty())
        std::cout << std::endl;
    if (json_out)
        std::cout << "}" << std::endl; // close everything
    return EXIT_SUCCESS;
}
