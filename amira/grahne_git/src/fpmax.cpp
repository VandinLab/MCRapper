/*
   Author:  Jianfei Zhu
			Concordia University
   Date:	Sep. 26, 2003

Copyright (c) 2003, Concordia University, Montreal, Canada
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   - Redistributions of source code must retain the above copyright notice,
	   this list of conditions and the following disclaimer.
   - Redistributions in binary form must reproduce the above copyright
	   notice, this list of conditions and the following disclaimer in the
	   documentation and/or other materials provided with the distribution.
   - Neither the name of Concordia University nor the names of its
	   contributors may be used to endorse or promote products derived from
	   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

*/

/* This is an implementation of FP-growth* / FPmax* /FPclose algorithm.
 *
 * last updated Sep. 26, 2003
 *
 */

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

#include <unistd.h>

#include "fsout.h"
#include "main.h"

void usage(const char *binary) {
    std::cerr << binary << ": mine frequent (closed|maximal) itemsets wrt "
        "minsupp from dataset" << std::endl
        << "USAGE: " << binary << " [-fhj] minsupp dataset" << std::endl
        << "\t-f : print additional information at the end of execution"
        << std::endl
        << "\t-h : print this message and exit" << std::endl
        << "\t-j : JSON output" << std::endl;
}

int main(int argc, char **argv) {
    bool full {false};
    bool json {false};
    {
        char opt;
        while ((opt = getopt(argc, argv, "fjh")) != -1) {
            switch (opt) {
                case 'f':
                    full = true;
                    break;
                case 'h':
                    usage(argv[0]);
                    return EXIT_SUCCESS;
                    break;
                case 'j':
                    json = true;
                    break;
            }
        }
    }
    if (optind != argc - 2) {
        std::cerr << "Error: wrong number of arguments" << std::endl;
        return EXIT_FAILURE;
    }
    const auto minsupp {argv[argc - 2]};
    const auto ds {argv[argc - 1]};

    const auto start {std::chrono::system_clock::now()};
    int rv;
    if (json) {
        std::cout << "{" << std::endl;
        std::string tabs {"\t"};
        if (full) {
            tabs = "\t\t";
            std::cout << "run: {" << std::endl
                << "\titemsets: {" << std::endl;
        }
        auto ftor {fpgrowth::JFSout(std::move(tabs))};
        rv = do_main(ds, minsupp, ftor);
    } else  {
        if (full) {
            std::cout << "# Run" << std::endl
                << std::endl << "## Itemsets" << std::endl;
        }
        auto ftor {fpgrowth::FSout()};
        rv = do_main(ds, minsupp, ftor);
    }
    const auto end {std::chrono::system_clock::now()};
    if (rv != EXIT_SUCCESS)
        return rv;
    if (full) {
        std::string comma;
        std::string quotes;
        std::string tab;
        if (json) {
            comma = ",";
            quotes = "\"";
            tab = "\t";
            std::cout << std::endl << "\t}," << std::endl // close itemsets
            << tab << "runtimes: {" << std::endl;
        } else
            std::cout << std::endl << "## Runtimes (ms)" << std::endl;
        std::cout << tab << tab << "total: " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    end - start).count() << std::endl;
        if (json)
            std::cout << tab << "}" << std::endl // close runtimes
                << "}," << std::endl // close run
                << "settings: {" << std::endl;
        else
            std::cout << std::endl << "# Settings" << std::endl;
        std::cout << tab << "algorithm: " << quotes << "FPGROWTH" << quotes
            << comma << std::endl
            << tab << "dataset: " << quotes << ds << quotes << comma
            << std::endl
            << tab << "minsupp: " << minsupp << std::endl;
        if (json) {
            std::cout << "}" << std::endl // close settings
                << "}" << std::endl; // close everything
        }
    } else if (json)
        // to close everything. The endline is needed because we are closing the
        // itemsets, and the last itemsets is printed with no newline after it.
        std::cout << std::endl << "}" << std::endl;;
    return EXIT_SUCCESS;
}
