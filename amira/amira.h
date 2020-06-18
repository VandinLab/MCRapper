/**
 * Useful aliases, functions, and classes for the AMIRA code.
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

#ifndef AMIRA_H_
#define AMIRA_H_

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <stdexcept>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

//#include <boost/container_hash/hash.hpp>
#include <boost/functional/hash.hpp>

namespace amira {

using count = unsigned long; // for sizes of sets of transactions (e.g., supports, sample sizes)
using item = unsigned int;
using itemset = std::vector<item>;

constexpr int meta {'#'};
//constexpr auto isuppmeta {"# isupp:"};
//constexpr auto sizemeta {"# size:"};
//constexpr auto isuppmetalen {8};
//constexpr auto sizemetalen {7};
constexpr std::string_view isuppmeta {"# isupp:"};
constexpr std::string_view sizemeta {"# size:"};
constexpr auto isuppmetalen {8};
constexpr auto sizemetalen {7};

// Compute the number of transactions in a file, using the metadata info if
// requested and found.
count get_size(const std::string &file, const bool use_metadata = true) {
    count size {0};
    std::ifstream in {file};
    if (in.fail())
        throw std::runtime_error(std::strerror(errno));
    for (std::string line; std::getline(in, line); ) {
        if (line[0] == meta) { // Is this line a metadata line?
            if (use_metadata && line.find(sizemeta) == 0) {
                size = std::stoul(line.substr(sizemetalen));
                break;
            }
        } else
            ++size;
    }
    return size;
}

// The following class stores the information we need for each itemset. The
// names of the members, if not documented, are taken from the paper.
struct ItemsetInfo {
    ItemsetInfo(const count sp) : sp {sp} {}

    // TODO: The use of mutable (and having update() being const) is ugly and it
    // is due to a wrong design. We should have classes for itemset,
    // ItemsetWithSupport, and ItemsetInfo (which should really be AmiraInfo).
    // This design would allow us to have maps that have keys of type
    // ItemsetWithSupport and values ItemsetInfos (i.e., AmiraInfos), and avoid
    // having these mutable and const attributes.
    // TODO: Do we need an ordered map for g?
    // // g has no key '0' because this information is useless.
    mutable std::map<std::size_t,count> g;
    // h[i] stores the value for h_{A,i}, with the exception of h[0], which is a
    // non-informative value which is never used.
    mutable std::vector<count> h;
    // w[i] stores the value for w_{A,i}, with the exception of w[0], which is a
    // non-informative value which is never used.
    mutable std::vector<count> w;
    count sp; // support

    void update(const count k, const count c) const {
        if (k > 0) {
            if (const auto g_it {g.lower_bound(k)};
                    g_it == g.cend() || g_it->first != k)
                g.emplace_hint(g_it, k, 1);
            else
                ++(g_it->second);
            if (k < h.size())
                // XXX MR: A man can dream of compilers that support C++17
                // features
                //std::for_each_n(h.begin() + 1, k + 1, [](int &x){ ++x; });
                std::for_each(h.begin() + 1, h.begin() + k + 1,
                        [](amira::count &x) { ++x; });
            else {
                std::for_each(h.begin(), h.end(),
                        [](amira::count &x) { ++x; });
                h.resize(k + 1, 1);
            }
            if (c > 0) {
                if (k < w.size())
                    // XXX MR: A man can dream of compilers that support C++17
                    // features
                    //std::for_each_n(w.begin() + 1, k + 1,
                    //                  [c](int &x){ x += c; });
                    std::for_each(w.begin() + 1, w.begin() + k + 1,
                            [c](amira::count &x) { x += c; });
                else {
                    std::for_each(w.begin(), w.end(),
                            [c](amira::count &x) { x += c; });
                    w.resize(k + 1, c);
                }
            }
        }
    }
};

// Helper class to associate an itemset with its info. We use it to store the
// CFIs.
class ItemsetWithInfo {
public:
    ItemsetWithInfo(const itemset &itms, const ItemsetInfo &info) : itms {itms},
        info {info} {}

    const itemset itms;
    ItemsetInfo info;
};

// Helper class to sort the items in an itemset according to their support in
// the sample.
// XXX MR: is it used anywhere?
class ItemSuppComp {
private:
    const std::vector<count> &itemsupps;
public:
    ItemSuppComp(const std::vector<count> &itemsupps) : itemsupps {itemsupps} {}

    bool operator()(const item a, const item b) const {
        return itemsupps[a] < itemsupps[b];
    }
};

// Functor to compute the "lesser" of two ItemsetsWithInfo, first by length of
// the itemset (shorter is "less") then by support (smaller support is
// "less"), then alphabetically ('a' is less than 'b').
// Not the same as the order <_H.
class LengthThenSuppComp {
public:
    bool operator()(const ItemsetWithInfo &a, const ItemsetWithInfo &b) const {
        if (a.itms.size() != b.itms.size())
            return a.itms.size() < b.itms.size();
        if (a.info.sp != b.info.sp)
            return a.info.sp < b.info.sp;
        if (a.itms.size() != b.itms.size())
            return a.itms.size() > b.itms.size();
        auto a_it {a.itms.cbegin()};
        auto b_it {b.itms.cbegin()};
        for (; a_it != a.itms.cend(); ++a_it, ++b_it)
            if (*a_it != *b_it)
                return *a_it < *b_it;
        return false;
    }
};

// Functor to compute the "lesser" of two ItemsetsWithInfo, first by support
// (smaller support is "less") then in inverted order by length (longer length
// is "less"), then alphabetically ('a' is less than 'b').
// This is the order <_H.
class SuppThenInvByLengthComp {
public:
    bool operator()(const ItemsetWithInfo &a, const ItemsetWithInfo &b) const {
        if (a.info.sp != b.info.sp)
            return a.info.sp < b.info.sp;
        if (a.itms.size() != b.itms.size())
            return a.itms.size() > b.itms.size();
        auto a_it {a.itms.cbegin()};
        auto b_it {b.itms.cbegin()};
        for (; a_it != a.itms.cend(); ++a_it, ++b_it)
            if (*a_it != *b_it)
                return *a_it < *b_it;
        return false;
    }
};

// Return a string with the elements in [b,e) separated by sep.
template<typename It> std::string write_container(const It b, const It e, const
        std::string_view &sep = " ") {
    if (b != e) {
        std::ostringstream oss;
            // Convert all but the last item to avoid a trailing separator.
        std::copy(b, e - 1, std::ostream_iterator<
                typename std::iterator_traits<It>::value_type>(oss, sep.data()));
        // Now add the last item with no delimiter following it
        oss << *(e - 1);
        return oss.str();
    }
    return "";
}

// Write an itemset and its support to std::cout, using the appropriate
// separators and spacing.
class ItemsetPrinter {
public:
    ItemsetPrinter(const std::string_view &comma, const std::string_view &prespace,
            const std::string_view &quotes, const std::string_view &sep,
            const std::string_view &supp_close, const std::string_view &supp_open) :
        comma {comma}, prespace {prespace}, quotes {quotes}, sep {sep},
              supp_close {supp_close}, supp_open {supp_open} {}

    void operator()(const amira::itemset &i, const amira::count sp) {
        if (printed_++ > 0)
            std::cout << comma << std::endl;
        std::cout << prespace << quotes
            << write_container(i.cbegin(), i.cend(), sep) << quotes << supp_open
            << sp << supp_close;
    }

    count printed() const { return printed_; }
private:
    count printed_ {0};
    const std::string_view &comma;
    const std::string_view &prespace;
    const std::string_view &quotes;
    const std::string_view &sep;
    const std::string_view &supp_close;
    const std::string_view &supp_open;
};


// Compute the hash of the itemset.
class ItemsetHash {
public:
    std::size_t operator()(const itemset &i) const { return boost::hash_value(i); }
};

// Obtain the FIs from the CFIs in the container cfis (which is a set of
// ItemsetWithInfo, but possibly sorted weirdly). Apply the functor ftor to
// each FI. The Ftor should take a const itemset & and a const count as
// arguments.
// The algorithm to obtain the CFIs from the FIs is Algorithm 6.4 of Tan,
// Steinback, Kumar, Introduction to Data Mining.
template<typename S, typename Ftor> void cfis_to_fis(const S &cfis,
        Ftor &ftor) {
    // Sort the cfis by length and then by supp, which is more appropriate for
    // this algorithm. We traverse the set in _reverse_ order, so we first visit
    // the longest cfis, and do so in decreasing order of support, and then go
    // to shorter cfis, and visit them in decreasing order of support, and so
    // on.
    std::set<ItemsetWithInfo, LengthThenSuppComp> cfiss {cfis.cbegin(),
        cfis.cend()};
    std::vector<std::pair<itemset,count>> f;
    auto it {cfiss.rbegin()};
    // Find all frequent itemsets of maximum size.
    auto k {it->itms.size()};
    do {
        // The sorting (also in later lines) is important for computing the
        // hashes for the unordered_set seen.
        itemset cp {it->itms.cbegin(), it->itms.cend()};
        std::sort(cp.begin(), cp.end());
        f.emplace_back(cp, it->info.sp);
        ftor(cp, it->info.sp);
    } while (++it != cfiss.rend() && it->itms.size() == k);
    do {
        --k;
        decltype(f) nf;
        std::unordered_set<itemset,ItemsetHash> seen;
        for (const auto &isp : f) {
            // Insert into seen and nf and call ftor for all the CFIs of size k
            // that have a support at least as large as the one of isp.first,
            // i.e., isp.second. CFIs of size k with support smaller than
            // isp.second cannot be subsets of isp.first, thus won't be
            // generated in the subsequent loop.
            for (; it != cfiss.rend() && it->itms.size() == k &&
                    it->info.sp >= isp.second; ++it) {
                itemset cp {it->itms.cbegin(), it->itms.cend()};
                std::sort(cp.begin(), cp.end());
                seen.insert(cp);
                nf.emplace_back(cp, it->info.sp);
                ftor(cp, it->info.sp);
            }
            itemset i {isp.first.cbegin() + 1, isp.first.cend()};
            item r {isp.first[0]};
            for (item idx {0}; idx <= i.size(); ++idx) {
                if (seen.insert(i).second) {
                    nf.emplace_back(i, isp.second);
                    ftor(i, isp.second);
                }
                if (idx < i.size())
                    std::swap(r, i[idx]);
            }
        }
        // Insert into nf and call ftor for all the CFIs of size k that have a
        // support smaller than all the itemsets in f.
        for (; it != cfiss.rend() && it->itms.size() == k; ++it) {
            itemset cp {it->itms.cbegin(), it->itms.cend()};
            std::sort(cp.begin(), cp.end());
            nf.emplace_back(cp, it->info.sp);
            ftor(cp, it->info.sp);
        }
        f = std::move(nf);
    } while (k > 1);
}

}  // namespace amira

#endif // AMIRA_H_
