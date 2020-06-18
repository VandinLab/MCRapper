/**
 * Functions to create, mine, and write a sample.
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

#ifndef _SAMPLE_H_
#define _SAMPLE_H_

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "wheelwright/matteoutils/Sampler.h"

#include "grahne/common.h"
#include "grahne/algos.h"

#include "amira.h"

namespace amira {

// Sample sample_size transactions with replacement from dataset, which has size
// ds_size. Store in ``sample'' the unique sampled transactions together with
// the number of times they appear in the sample. The sum of the values of
// ``sample'' is sample_size.  The information infos needed for computing the
// first error bounds are also populated.
void create_sample(const std::string &dataset, const count ds_size,
        const count sample_size,
        std::unordered_map<itemset, count, ItemsetHash> &sample,
        std::map<item,ItemsetInfo> &infos) {
    {
        std::ifstream ds {dataset};
        if (ds.fail())
            throw std::runtime_error(std::strerror(errno));
        // indices (in the dataset) of the transactions to sample
        std::vector<std::size_t> tosample_idxs;
        if (sample_size != ds_size) {
            tosample_idxs = MatteoUtils::Sampler(
                static_cast<std::size_t>(ds_size)).samples(sample_size);
            std::sort(tosample_idxs.begin(), tosample_idxs.end());
        } else { // Consider the dataset as the sample.
            tosample_idxs.resize(ds_size);
            std::iota(tosample_idxs.begin(), tosample_idxs.end(), 0);
        }
        // index in the dataset of the transaction just read
        std::size_t curr_idx {0};
        auto idxsit {tosample_idxs.begin()};
        for (std::string line; idxsit != tosample_idxs.cend() &&
                std::getline(ds, line);) {
            // Skip metadata lines.
            if (line.at(0) == amira::meta)
                continue;
            // Go to the next transaction if the transaction just read should
            // not be sampled.
            if (curr_idx != *idxsit) {
                ++curr_idx;
                continue;
            }
            std::istringstream ls {line};
            itemset t {std::istream_iterator<item>(ls),
                std::istream_iterator<item>()};
            // We sort the items in a transaction in alphabetical order because
            // we need a consistent ordering of items in transactions and cfis
            // to populate the itemset infos needed to compute omega2 and rho2.
            std::sort(t.begin(), t.end());
            count copies {0};
            // Add the transaction to the sample (possibly multiple times).
            do {
                ++idxsit;
                ++copies;
            } while (idxsit != tosample_idxs.end() && curr_idx == *idxsit);
            for (const auto i : t) {
                if (const auto infit {infos.lower_bound(i)};
                        infit == infos.end() || infit->first != i)
                    infos.emplace_hint(infit, i, copies);
                else
                    infit->second.sp += copies;
            }
            if (auto it {sample.find(t)}; it == sample.cend())
                sample.emplace(std::move(t), std::move(copies));
            else
                it->second += copies;
            ++curr_idx;
        }
    }
    // Populate the infos of the items appearing in the sample.
    // The members g and h of the ItemsetInfo's in infos are only updated the
    // for each distinct transaction, because copies do not change the number of
    // closed itemsets.
    for (const auto &p : sample) {
        // Find the item in t that comes first wrt the order <_H (i.e., the
        // item with minimum support, and in case of equal support, the
        // smallest item).
        const auto m {std::min_element(p.first.cbegin(), p.first.cend(),
                [&infos = std::as_const(infos)](const item a, const item b) {
                    if (infos.at(a).sp == infos.at(b).sp)
                        return a < b;
                    else
                        return infos.at(a).sp < infos.at(b).sp;
                })};
        infos.at(*m).update(p.first.size() - 1, p.second);
    }
}

// Mine the closed frequent itemsets in the sample s at the minimum support
// threshold thres, and call ftor on each of them.
template<typename Ftor> void mine_sample(
        const std::unordered_map<itemset, count, ItemsetHash> &s,
        const count thres, const std::map<item, ItemsetInfo> &infos,
        Ftor &ftor) {
    // In order to use the FPGrowth code by Grahne and Zhu without too many
    // modifications, we need to collect the item supports in a vector
    // supps such that supps[i] is the support of item i.
    std::vector<int> supps(infos.rbegin()->first + 1);
    for (const auto &it : infos)
        supps[it.first] = it.second.sp;
    // Most of the following code is taken from grahne/main.h .
    fpgrowth::THRESHOLD = thres;
    fpgrowth::FI_tree* fptree;
    fpgrowth::fp_buf = new fpgrowth::memory(60, 4194304L, 8388608L, 2);
    fptree = (fpgrowth::FI_tree*)fpgrowth::fp_buf->newbuf(1,
            sizeof(fpgrowth::FI_tree));
    fptree->init(-1, 0);
    fptree->scan1_DB_init_arrays(supps, infos.crbegin()->first);
    fpgrowth::ITlen = new int[fptree->itemno] {}; // C++11 default init
    fpgrowth::bran = new int[fptree->itemno] {}; // C++11 default init
    fpgrowth::compact = new int[fptree->itemno];
    fpgrowth::prefix = new int[fptree->itemno];
    fpgrowth::list = new fpgrowth::stack(fptree->itemno, true);
    fptree->scan2_DB(s);
    // We skip handling the emptyset, because we don't need it
    if(fptree->Single_path()) {
        int i {0};
        for(fpgrowth::Fnode *node {fptree->Root->leftchild}; node != NULL;
                node = node->leftchild) {
            fpgrowth::list->FS[i] = node->itemname;
            fpgrowth::list->counts[i++] = node->count;
        }
        int Count;
        i = 0;
        while (i < fptree->itemno) {
            Count = fpgrowth::list->counts[i];
            for(; i < fptree->itemno && fpgrowth::list->counts[i] == Count;
                    i++)
                ; /* empty for loop */
            // i is necessarily at least 1, because the above for loop is
            // executed at least once.
            fpgrowth::ITlen[i-1]++;
            ftor(i, fpgrowth::list->FS, Count);
        }
        return;
    }
    fpgrowth::current_fi = new bool[fptree->itemno] {}; // C++11 default init
    fpgrowth::supp = new int[fptree->itemno] {}; //for keeping support of items
    fpgrowth::CFI_tree* LClose;
    fpgrowth::cfitrees =
        (fpgrowth::CFI_tree**) new fpgrowth::CFI_tree*[fptree->itemno];
    fpgrowth::memory* Close_buf =
        new fpgrowth::memory(40, 1048576L, 5242880, 2);
    LClose =
        (fpgrowth::CFI_tree*) Close_buf->newbuf(1, sizeof(fpgrowth::CFI_tree));
    LClose->init(Close_buf, fptree, NULL, NULL, -1);
    fptree->set_close_tree(LClose);
    fpgrowth::cfitrees[0] = LClose;
    fptree->FPclose(ftor);

    delete Close_buf;
    delete []fpgrowth::cfitrees;
    delete []fpgrowth::supp;
    delete []fpgrowth::current_fi;
    delete fpgrowth::list;
    delete []fpgrowth::prefix;
    delete []fpgrowth::compact;
    delete []fpgrowth::bran;
    delete []fpgrowth::ITlen;
    delete fpgrowth::fp_buf;
    delete []fpgrowth::order_item;
    delete []fpgrowth::item_order;
}

void write_sample(const std::unordered_map<itemset, count, ItemsetHash> &s,
        const std::string &out) {
    std::ofstream of {out};
    if (of.fail())
        throw std::runtime_error(std::strerror(errno));
    for (const auto &p : s) {
        const auto w {amira::write_container(p.first.cbegin(), p.first.cend())};
        for (count i {0}; i < p.second; ++i)
            of << w << std::endl;
    }
}

// Functor used to get the results of the mining and store them in the set cfis,
// sorted according to the ordering <_H.
// This class is defined here because it accesses structures in the fpgrowth
// namespace, and we don't want pollution other files with headers too much.
class AddItemsetToSet {
    public:
        AddItemsetToSet(std::set<ItemsetWithInfo,
                SuppThenInvByLengthComp> &cfis) : cfis {cfis} {}
        void operator()(const int ilen, const int *iset, const int support,
          const int plen = 0, const int *prefix = nullptr) {
            if (int sz {ilen + plen}; sz > 0) {
                amira::itemset itms;
                itms.reserve(sz);
                for (int i {0}; i < plen; ++i)
                    itms.push_back(fpgrowth::order_item[prefix[i]]);
                for (int i {0}; i < ilen; ++i)
                    itms.push_back(fpgrowth::order_item[iset[i]]);
                // We sort the items in the itemset in alphabetical order
                // because we need a consistent ordering of itemset in
                // transactions and cfis to populate the itemset infos needed to
                // compute omega2 and rho2.
                std::sort(itms.begin(), itms.end());
                cfis.emplace(itms, static_cast<amira::count>(support));
            }
        }

    private:
        std::set<ItemsetWithInfo, SuppThenInvByLengthComp> &cfis;
};


} // namespace amira
#endif // _SAMPLE_H_
