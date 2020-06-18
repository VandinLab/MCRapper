#ifndef _FSOUT_H_
#define _FSOUT_H_

#include <iostream>
#include <string>

#include "common.h"

namespace fpgrowth {

class FSout {
public:
    void operator()(const int ilen, const int *iset, const int support,
            const int plen = 0, const int *prefix = nullptr) {
        for(int i {0}; i < plen; ++i)
            std::cout << order_item[prefix[i]] << " ";
        for(int i {0}; i < ilen; ++i)
            std::cout << order_item[iset[i]] << " ";
        std::cout << "(" << support << ")" << std::endl;
    }
};

class JFSout {
public:
    JFSout(const std::string &tabs = "") : tabs {tabs} {}
    JFSout(std::string &&tabs) : tabs {tabs} {}

    void operator()(const int ilen, const int *iset, const int support,
            const int plen = 0, const int *prefix = nullptr) {
        if (notfirst)
            std::cout << "," << std::endl;
        else
            notfirst = true;
        std::cout << tabs << "\"";
        if (ilen + plen == 0)
            std::cout << "*";
        else {
            for(int i {0}; i < plen; ++i)
                std::cout << order_item[prefix[i]] << "_";
            for(int i {0}; i < ilen-1; ++i)
                std::cout << order_item[iset[i]] << "_";
            if (ilen > 0)
                std::cout << order_item[iset[ilen-1]];
        }
        std::cout << "\": " << support;
    }
private:
    bool notfirst {false};
    const std::string tabs;
};


class DoNothingFtor {
public:
    void operator()(const int, const int *, const int,
            const int = 0, const int* = nullptr) {}
};

} // namespace fpgrowth

#endif // _FSOUT_H_
