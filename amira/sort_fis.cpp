/**
 * Sort a collection of itemsets in decreasing order according to their
 * support, and sort the item in each itemset according to their numerical
 * value. The input to this utility must be the output of either amira or
 * grahne/fim_all (run with or without the -f and/or -j flags).
 *
 * Copyright 2014, 2018 Matteo Riondato <riondato@acm.org>
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

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>

#include "amira.h"
#include "amira_regex.h"

void usage(const char *binary) {
    std::cerr << binary << ": sort itemsets in decreasing order of support"
        << std::endl
        << "USAGE: " << binary << " [-hjv] fis" << std::endl
        << "\t-h : print this message and exit" << std::endl
        << "\t-j : JSON input (i.e., amira or fim_* were run with the "
        "'-j' flag)" << std::endl
        << "\t-v : print log messages to stderr during the execution"
        << std::endl;
}

std::string sort_itms(const std::string &str, const std::string_view &sep) {
    auto start {str.c_str()};
    char *end;
    std::vector<amira::item> itms;
    auto itm {std::strtoul(start, &end, 10)};
    for (; start != end; itm = std::strtoul(start, &end, 10)) {
        start = end;
        itms.emplace_back(itm);
    }
    std::sort(itms.begin(), itms.end());
    return amira::write_container(itms.cbegin(), itms.cend(), sep);
}

int main(int argc, char **argv) {
    bool json {false};
    bool verbose {false};
    {
        extern int optind;
        char opt;
        while ((opt = getopt(argc, argv, "hjv")) != -1) {
            switch (opt) {
                case 'h':
                    usage(argv[0]);
                    return EXIT_SUCCESS;
                    break;
                case 'j':
                    json = true;
                    break;
                case 'v':
                    verbose = true;
                    break;
                default:
                    std::cerr << "Error: wrong option." << std::endl;
                    return EXIT_FAILURE;
            }
        }
        if (optind != argc - 1) {
            std::cerr << "Error: wrong number of arguments" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::ifstream f {argv[argc - 1]};
    if (f.fail()) {
        std::cerr << "Error: can't open the fis file." << std::endl;
        return EXIT_FAILURE;
    }
    std::vector<std::pair<std::string,amira::count>> s;
    std::string line;
    std::cmatch m;
    std::string_view around_itms;
    std::string_view open_supp {"("};
    std::string_view close_supp {")"};
    std::string_view comma;
    std::string_view sep {" "};
    if (json) {
        around_itms = "\"";
        open_supp = ": ";
        close_supp = "";
        comma = ",";
        sep = "_";
        if (verbose)
            std::cerr << "Looking for itemsets in JSON input...";
        // Print the heading invalid lines
        while (std::getline(f, line) &&
                !  std::regex_match(line.c_str(), m, amira::JSONFIre))
            std::cout << line << std::endl;
        if (! f.fail()) {
            // line contains the first valid line. Process it and continue
            // processing as long as we get valid lines.
            do {
                std::string itms {m[1].first, m[1].second};
                if (itms != "*") {
                    std::replace(itms.begin(), itms.end(), '_', ' ');
                    itms = sort_itms(itms, sep);
                }
                s.emplace_back(itms, std::strtoul(m[2].first, NULL, 10));
            } while (std::getline(f, line) &&
                    std::regex_match(line.c_str(), m, amira::JSONFIre));
        }
    } else {
        if (verbose)
            std::cerr << "Looking for itemsets in plaintext input...";
        // Print the heading invalid lines
        while (std::getline(f, line) &&
                !  std::regex_match(line.c_str(), m, amira::FIre))
            std::cout << line << std::endl;
        if (! f.fail()) {
            // line contains the first valid line. Process it and continue
            // processing as long as we get valid lines.
            do {
                std::string itms {m[1].first, m[1].second};
                if (! itms.empty())
                    itms = sort_itms(itms, sep) + " ";
                s.emplace_back(itms, std::strtoul(m[2].first, NULL, 10));
            } while (std::getline(f, line) &&
                    std::regex_match(line.c_str(), m, amira::FIre));
        }
    }
    if (verbose)
        std::cerr << "done: " << s.size() << " itemsets" << std::endl;
    if (! s.empty()) {
        if (verbose)
            std::cerr << "Sorting...";
        // Sort the itemsets in decreasing order of support, breaking ties by
        // increasing length of the itemsets, with the exception of the empty
        // itemset in the JSON case, which must be the first in the order.
        std::sort(s.rbegin(), s.rend(),
                [json](const std::pair<std::string,amira::count> &a,
                    const std::pair<std::string,amira::count> &b) -> auto {
                if (a.second != b.second)
                    return a.second < b.second;
                else if (json && a.first.size() == 1 && b.first.size() == 1 &&
                        (a.first == "*" || b.first == "*")) {
                    return a.first != "*";
                } else
                    return ! (a.first.size() < b.first.size());
                });
        if (verbose)
            std::cerr << "done" << std::endl << "Printing output and exiting. "
                "Goodbye." << std::endl;
        // Print the sorted itemsets
        for (auto it {s.cbegin()}; it != s.cend() - 1; ++it)
            std::cout << around_itms << it->first << around_itms << open_supp
                << it->second << close_supp << comma << std::endl;
        std::cout << around_itms << s.rbegin()->first << around_itms
            << open_supp << s.rbegin()->second << close_supp << std::endl;
    }
    // Print the remaining invalid lines;
    if (! f.fail()) {
        do
            std::cout << line << std::endl;
        while (std::getline(f, line));
    }
    return EXIT_SUCCESS;
}
