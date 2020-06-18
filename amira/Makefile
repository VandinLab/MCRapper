flags=-O2 -ftree-loop-vectorize -march=native -mtune=native -pipe -std=gnu++17 -fopenmp
warnsflags=-Wall -Wextra -Wdouble-promotion -Wformat=2 -Wnull-dereference -Wswitch-default -Wswitch-enum -Wswitch-bool -Wswitch-unreachable -Wsync-nand -Wuninitialized -Wsuggest-final-types -Wsuggest-final-methods -Wsuggest-override -Wduplicated-branches -Wduplicated-cond -Wundef -Wcast-qual -Wlogical-op -Wextra-semi -Wredundant-decls -Wvirtual-move-assign -pedantic
debugflags=-ggdb3 -grecord-gcc-switches -D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables
linkingflags=-lm -lnlopt

CPP=g++
os=$(shell sh -c 'uname 2>/dev/null || echo Unknown')
ifeq ($(os), Darwin)
CPP=g++-mp-9
else
ifeq ($(os), FreeBSD)
CPP=g++8
ltoflags=-flto
addflags=-Wl,-rpath=/usr/local/lib/gcc8
else
ifeq ($(os), Linux)
CPP=g++-8
ltoflags=-flto
endif
endif
endif

all: amira compare_fis sort_fis

amira: Makefile amira.cpp amira.h epsilon.h sample.h grahne/algos.h grahne/buffer.cpp grahne/buffer.h grahne/common.h grahne/fp_node.cpp grahne/fp_node.h grahne/fp_tree.cpp grahne/fp_tree.h grahne/glob.cpp
	$(CPP) -o amira -DWHEELWRIGHT_NLOPT $(flags) $(ltoflags) $(warnsflags) $(debugflags) amira.cpp grahne/buffer.cpp grahne/fp_node.cpp grahne/fp_tree.cpp grahne/glob.cpp $(linkingflags) $(addflags)

compare_fis: Makefile compare_fis.cpp amira.h amira_regex.h
	$(CPP) -o compare_fis $(flags) $(ltoflags) $(warnsflags) $(debugflags) compare_fis.cpp $(linkingflags) $(addflags)

sort_fis: Makefile sort_fis.cpp amira.h amira_regex.h
	$(CPP) -o sort_fis $(flags) $(warnsflags) $(ltoflags) $(debugflags) sort_fis.cpp $(linkingflags) $(addflags)

.PHONY: clean
clean:
	-rm -rf amira compare_fis sort_fis *.dSYM
