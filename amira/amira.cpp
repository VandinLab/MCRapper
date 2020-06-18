/**
 * The AMIRA algorithm for fixed-size samples.
 *
 * Copyright 2018-2019 Matteo Riondato <riondato@acm.org>
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
#include <cerrno>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <unistd.h>

#include "amira.h"
#include "epsilon.h"
#include "sample.h"

#include "grahne/common.h"

void usage(const char *binary, const int code) {
    std::cerr << binary << ": run AMIRA with a fixed sample size" << std::endl
        << "USAGE: " << binary << " [-c] [-d dataset_size] [-fhjnp] "
        "[-r maximum_frequency] [-s sample] [-v] failure_probability "
        "minimum_frequency sample_size dataset" << std::endl
        << "\t-c : print the closed frequent itemsets, rather than the "
        "frequent ones" << std::endl
        << "\t-d : specify the dataset size as 'dataset_size'" << std::endl
        << "\t-f : print full information about the run at the end" << std::endl
        << "\t-h : print this message and exit" << std::endl
        << "\t-j : print final output in JSON format" << std::endl
        << "\t-n : do not output the itemsets at the end" << std::endl
        << "\t-p : skip the computation of the second bound" << std::endl
        << "\t-r : assume items have frequency at most 'maximum_frequency'"
        << std::endl
        << "\t-s : write the sampled transactions to file 'sample'" << std::endl
        << "\t-v : print log messages to stderr during the execution"
        << std::endl;
    std::exit(code);
}

int main(int argc, char** argv) {
    bool full {false};
    bool json {false};
    bool noitmsets {false};
    bool printclosed {false};
    bool skipsecond {false};
    bool verbose {false};
    double maxfreq {1};
    amira::count ds_size {0};
    std::string outf;
    char opt;
    extern char *optarg;
    extern int optind;
    while ((opt = getopt(argc, argv, "cd:fhjnpr:s:v")) != -1) {
        switch (opt) {
            case 'c':
                printclosed = true;
                break;
            case 'd':
                ds_size = std::strtoul(optarg, NULL, 10);
                if (errno == ERANGE || errno == EINVAL || ds_size == 0) {
                    std::cerr << "Error: ds_size must be a positive integer"
                        << std::endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'f':
                full = true;
                break;
            case 'h':
                usage(argv[0], EXIT_SUCCESS);
                break;
            case 'j':
                json = true;
                break;
            case 'n':
                noitmsets = true;
                break;
            case 'p':
                skipsecond = true;
                break;
            case 'r':
                maxfreq = std::strtod(optarg, NULL);
                if (errno == ERANGE || maxfreq <= 0 || maxfreq > 1) {
                    std::cerr << "Error: maximum_frequency must be a real in "
                        "(0,1]" << std::endl;
                    return EXIT_FAILURE;
                }
                break;
            case 's':
                outf = std::string(optarg);
                break;
            case 'v':
                verbose = true;
                break;
            default:
                std::cerr << "Error: wrong option." << std::endl;
                return EXIT_FAILURE;
        }
    }
    if (optind != argc - 4) {
        std::cerr << "Error: wrong number of arguments" << std::endl;
        return EXIT_FAILURE;
    }
    // Error probability
    const double delta {std::strtod(argv[argc - 4], NULL)};
    if (errno == ERANGE || delta <= 0 || delta >= 1) {
        std::cerr << "Error: delta must be a real in (0,1)" << std::endl;
        return EXIT_FAILURE;
    }
    // Original frequency threshold
    const double theta {std::strtod(argv[argc - 3], NULL)};
    if (errno == ERANGE || theta <= 0 || theta >= 1) {
        std::cerr << "Error: theta must be a real in (0,1)" << std::endl;
        return EXIT_FAILURE;
    }
    amira::count size {std::strtoul(argv[argc - 2], NULL, 10)};
    if (errno == ERANGE) {
        std::cerr << "Error: samplesize must be a positive integer"
            << std::endl;
        return EXIT_FAILURE;
    }
    // XXX MR: uncomment the following and remove the above conversion once
    // there is support for std::from_chars in compilers.
    //{
        // We use from_chars because it doesn't require us to think about what
        // type is 'size', so we can change the type without having to change
        // the code.
        //const auto last {argv[argc - 2] + std::strlen(argv[argc - 2])};
        //const auto r {std::from_chars(argv[argc - 2], last, size)};
        //if (r.ptr != last || size <= 0) {
        //    std::cerr << "Error: samplesize must be a positive integer"
        //        << std::endl;
        //    return EXIT_FAILURE;
        //}
    //}
    const std::string dataset {argv[argc - 1]};
    const auto start {std::chrono::system_clock::now()};
    if (ds_size == 0) {
        if (verbose)
            std::cerr << "Getting dataset size...";
        try {
            ds_size = amira::get_size(dataset);
        } catch (std::runtime_error &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        if (verbose)
            std::cerr << "done (" << ds_size << " transactions)" << std::endl;
    }
    if (verbose)
        std::cerr << "Creating sample of size " << size << "...";
    // The unique sampled transactions, with the number of times they appear in
    // the sample.
    std::unordered_map<amira::itemset, amira::count, amira::ItemsetHash> sample;
    std::map<amira::item,amira::ItemsetInfo> item_infos;
    // Create the sample and populate the item_infos data structure needed to
    // compute the first omega and the first rho.
    // XXX MR: Strictly speaking, in terms of running time, populating the data
    // structure should probably be assigned to the time needed to compute
    // omega1 and rho1, but implementation-wise, it is easier to populate the
    // structure in create_sample.
    try {
        // TODO MR: if notitmsets && skipsecond, we skip the mining completely,
        // so we can do _less_ in create_sample, so we should consider
        // implementing a special function create__sample_no_mining.
        amira::create_sample(dataset, ds_size, size, sample, item_infos);
    } catch (std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    // Set the support of items to be at most maxsupp
    const amira::count maxsupp = (amira::count) std::floor(size * maxfreq);
    for (auto &infit : item_infos) {
        if (infit.second.sp > maxsupp)
            infit.second.sp = maxsupp;
    }
    const auto create_sample_end {std::chrono::system_clock::now()};
    // Write the sample to an output file if requested.
    if (! outf.empty()) {
        if (verbose)
            std::cerr << "done" << std::endl << "Writing sample...";
        try {
            amira::write_sample(sample, outf);
        } catch (std::runtime_error &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
    }
    if (verbose)
        std::cerr << "done" << std::endl << "Computing omega1 and rho1...";
    // Compute omega1 and rho1 using the item supports.
    const auto rho1_start {std::chrono::system_clock::now()};
    amira::EraEps item_er; // stores omega1 and rho1
    try {
        amira::compute_eraeps(delta, size, item_infos, item_er);
    } catch (const std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    const auto rho1_end {std::chrono::system_clock::now()};
    if (verbose)
        std::cerr << "done (omega1: " << item_er.era << ", rho1: "
            << item_er.eps << ")" << std::endl;
    double freq1 {-1.0};
    amira::count supp1 {0};
    std::set<amira::ItemsetWithInfo, amira::SuppThenInvByLengthComp> q;
    auto mine_end {std::chrono::system_clock::now()};
    std::size_t cfis1 {0};
    amira::EraEps itemset_er {-1.0,-1.0}; // stores omega2 and rho2
    double r {item_er.eps};
    double freq2 {-1.0};
    amira::count supp2 {0};
    auto rho2_end {mine_end};
    auto prune_end {mine_end};
    if (noitmsets && skipsecond)
        goto output;
    { // limit the scope of variables only needed in case we perform the mining.
        if (verbose)
            std::cerr << "Mining CFIs from the sample at frequency ";
        // Compute the first lowered frequency and support thresholds.
        freq1 = theta - item_er.eps;
        if (freq1 <= 0)
            freq1 = 1.0 / size;
        if (verbose)
            std::cerr << freq1 << "...";
        supp1 = static_cast<amira::count>(std::ceil(freq1 * size));
        // Mine the Closed Frequent Itemesets (CFIs) in the sample at the first
        // lowered support threshold. The CFIs are stored in q, which is ordered
        // according to <_q.
        amira::AddItemsetToSet ftor {q};
        amira::mine_sample(sample, supp1, item_infos, ftor);
        mine_end = rho2_end = prune_end = std::chrono::system_clock::now();
        cfis1 = q.size();
        if (verbose)
            std::cerr << "done (" << cfis1  << " CFIs found)" << std::endl;
        amira::count supp {supp1};
        if (skipsecond)
            goto output;
        if (verbose)
            std::cerr << "Computing omega2 and rho2...";
        // Add to q the items that are not frequent wrt freq1
        {
            for (const auto &ii : item_infos)
                if (ii.second.sp < supp1)
                    // We cannot reuse ii.second because its g and h are already
                    // populated, while we want ``vanilla' ones.
                    q.emplace(amira::itemset({ii.first}),
                            amira::ItemsetInfo(ii.second.sp));
        }
        // Populate the infos of the itemsets in q: for each transaction t we
        // find the itemset in t that appears first in q.
        // XXX MR: An openMP for-loop would help here, but we need to figure out
        // the reduction or maybe the locks, as we write into q.
        for (const auto &p : sample) {
            for (auto qit {q.begin()}; qit != q.end(); ++qit) {
                // std::includes() requires consistent sorting between the
                // elements in the two ranges. We sorted the items in each
                // transaction at sample creation time and sorted the items
                // in each element of q in impl::AddItemsetToSet while
                // populating q.
                // TODO MR: Is std::includes slow? Shall we try rolling our
                // own using binary search?
                if (const auto& cfi {*qit}; cfi.itms.size() <= p.first.size() &&
                        std::includes(p.first.cbegin(), p.first.cend(),
                            cfi.itms.cbegin(), cfi.itms.cend())) {
                    qit->info.update(p.first.size() - qit->itms.size(),
                            p.second);
                    break;
                }
            }
        }
        // Compute omega2 and rho2
        try {
            const amira::count supp1_minus1 {supp1 - 1};
            auto data {std::tie(supp1_minus1, std::as_const(q))};
            amira::compute_eraeps(delta, size, data, itemset_er);
        } catch (const std::runtime_error &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        // Compute the second lowered frequency and support thresholds.
        freq2 = theta - itemset_er.eps;
        if (freq2 <= 0)
            freq2 = 1.0 / size;
        supp2 = static_cast<amira::count>(std::ceil(freq2 * size));
        if (itemset_er.eps < r) {
            r = itemset_er.eps;
            supp = supp2;
        }
        rho2_end = std::chrono::system_clock::now();
        if (verbose)
            std::cerr << "done (omega2: " << itemset_er.era << ", rho2: "
                << itemset_er.eps << ")" << std::endl;
        if (! noitmsets) {
            if (verbose)
                std::cerr << "Pruning CFIs...";
            // Remove from q the itemsets with support less than supp.
            q.erase(q.begin(), q.upper_bound({{}, supp - 1}));
        }
        prune_end = std::chrono::system_clock::now();
        if (verbose) {
            if (! noitmsets)
                std::cerr << "done (from " << cfis1 << " to " << q.size()
                    << " CFIs)" << std::endl;
        }
    } // end of reducing the scope of variables only needed for the mining.
output:
    if (verbose)
        std::cerr << "Printing output and exiting. Goodbye." << std::endl;
    // Printing output
    std::string_view comma;
    std::string_view quotes;
    std::string_view sep {" "};
    std::string_view supp_open {" ("};
    std::string_view supp_close {")"};
    std::string_view tab;
    std::string_view tabtab;
    if (json) {
        comma = ",";
        quotes = "\"";
        sep = "_";
        supp_close = "";
        supp_open = ": ";
        tab = "\t";
        tabtab = "\t\t";
        std::cout << "{" << std::endl; // open everything
    }
    if (full) {
        const auto starttime {std::chrono::system_clock::to_time_t(start)};
        std::cout << quotes << "date" << quotes << ": " << quotes
            << std::put_time(std::localtime(&starttime), "%F %T") << quotes
            << comma << std::endl;
        if (json)
            std::cout << quotes << "settings" << quotes << ": {" << std::endl;
        else
            std::cout << std::endl << "# Settings" << std::endl;
        std::cout << tab << quotes << "algorithm" << quotes << ": " << quotes
            << "AMIRA" << quotes << comma << std::endl
            << tab << quotes << "dataset" << quotes << ": " << quotes << dataset
            << quotes << comma << std::endl
            << tab << quotes << "samplesize" << quotes << ": " << size << comma
            << std::endl
            << tab << quotes << "minimum_frequency" << quotes << ": " << theta
            << comma << std::endl
            << tab << quotes << "failure_probability" << quotes << ": " << delta
            << comma << std::endl
            << tab << quotes << "maximum_frequency" << quotes << ": " << maxfreq
            << comma << std::endl
            << tab << quotes << "printclosed" << quotes << ": "
            << ((printclosed) ? 1 : 0) << comma << std::endl
            << tab << quotes << "skipsecond" << quotes << ": "
            << ((skipsecond) ? 1 : 0) << comma << std::endl
            << tab << quotes << "sample" << quotes << ": " << quotes
            << ((outf.empty()) ? "N/A" : outf) << quotes << std::endl;
        if (json) {
            // close settings
            std::cout << "}," << std::endl
                << quotes << "run" << quotes << ": {" << std::endl;
        } else
            std::cout << std::endl << "# Run" << std::endl;
        std::cout << tab << quotes << "omega1" << quotes << ": " << item_er.era
            << comma << std::endl;
        std::cout << tab << quotes << "rho1" << quotes << ": " << item_er.eps
            << comma << std::endl;
        std::cout << tab << quotes << "freq1" << quotes << ": " << freq1
            << comma << std::endl;
        std::cout << tab << quotes << "supp1" << quotes << ": " << supp1
            << comma << std::endl;
        std::cout << tab << quotes << "cfis1" << quotes << ": " << cfis1
            << comma << std::endl;
        std::cout << tab << quotes << "omega2" << quotes << ": "
            << itemset_er.era << comma << std::endl;
        std::cout << tab << quotes << "rho2" << quotes << ": " << itemset_er.eps
            << comma << std::endl;
        std::cout << tab << quotes << "freq2"<< quotes << ": " << freq2 << comma
            << std::endl;
        std::cout << tab << quotes << "supp2"<< quotes << ": " << supp2 << comma
            << std::endl;
        std::cout << tab << quotes << "cfis2"<< quotes << ": " << q.size()
            << comma << std::endl;
    } // full
    std::cout << tab << quotes << "eps" << quotes << ": " << 2*r;
    if (full || ! noitmsets)
        std::cout << comma;
    std:: cout << std::endl;
    if (full) {
        if (json)
            std::cout << tab << quotes << "runtimes" << quotes << ": {"
                << std::endl;
        else
            std::cout << std::endl << "## Runtimes (ms)" << std::endl;
        std::cout << tab << tab << quotes << "total" << quotes << ": " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    // We remove the time taken to write the sample to an
                    // output file.
                    (prune_end - start) - (rho1_start - create_sample_end)
                    ).count()
            << comma << std::endl
            << tab << tab << quotes << "create_sample" << quotes << ": " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    create_sample_end - start).count() << comma << std::endl
            << tab << tab << quotes << "get_rho1" << quotes << ": " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    rho1_end - rho1_start).count() << comma << std::endl
            << tab << tab << quotes << "mine" << quotes << ": " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    mine_end - rho1_end).count() << comma << std::endl
            << tab << tab << quotes << "get_rho2" << quotes << ": " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    rho2_end - mine_end).count() << comma << std::endl
            << tab << tab << quotes << "prune" << quotes << ": " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    prune_end - rho2_end).count() << std::endl;
        if (json)
            std::cout << tab << "}"; // close runtimes.
    } // full
    if (! noitmsets) {
        if (json) {
            if (full)
                std::cout << "," << std::endl;
            std::cout << tab << quotes << "itemsets" << quotes << ": {"
                << std::endl;
        } else
            std::cout << std::endl << "## Itemsets" << std::endl;
        std::string emptyset;
        std::string space;
        if (json) {
            emptyset = "*";
            space = " ";
        }
        amira::ItemsetPrinter ip {comma, tabtab, quotes, sep, supp_close,
            supp_open};
        if (printclosed) {
            bool printemptyset {true};
            for (const auto &e : item_infos) {
                if (e.second.sp == size) {
                    printemptyset = false;
                    break;
                }
            }
            if (printemptyset)
                std::cout << tab << tab << quotes << emptyset << quotes
                    << supp_open << size << supp_close << comma << std::endl;
            // Traverse q in reverse order to print cfis in decreasing order of
            // support.
            std::for_each(q.crbegin(), q.crend(), [&ip](const auto &cfi) {
                ip(cfi.itms, cfi.info.sp);});

        } else {
            std::cout << tab << tab << quotes << emptyset << quotes << supp_open
                << size << supp_close << comma << std::endl;
            // XXX MR: It is fair to ask whether the time taken in computing the
            // conversion from CFIs to FIs should be accounted for in our
            // runtime.We believe that it should not, because, even if only
            // intrinsically, the computation of the approximation is complete
            // by the time we have pruned q (which is not done if the user does
            // not ask for the itemsets, but that's a non-default variant).
            amira::cfis_to_fis(q, ip);
        }
        if (json) {
            std::cout << std::endl << "\t}"; // close itemsets
            if (full) {
                // Print the number of itemsets. Only printed in the json &&
                // full case, but it should be sufficient.
                std::cout << "," << std::endl << tab << quotes
                    << "itemsets_num" << quotes << ": " << ip.printed()
                    << std::endl << "}" << std::endl; // close run
            }
            std::cout << "}" << std::endl; // close everything
        } else
            std::cout << std::endl;
        // end of if( ! noitemsets)
    } else if (json) {
        // close run or everything (if !full)
        std::cout << std::endl << "}" << std::endl;
        if (full)
            std::cout << "}" << std::endl; // close everything
    }
    return EXIT_SUCCESS;
}
