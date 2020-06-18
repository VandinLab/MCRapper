#ifndef _MAIN_H_
#define _MAIN_H_

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "algos.h"
#include "common.h"
#include "buffer.h"
#include "fp_tree.h"

template<typename Ftor> int do_main(const char *dsf, const char *th, Ftor &ftor) {
    try {
        fpgrowth::THRESHOLD = std::stoi(std::string(th));
    } catch (const std::invalid_argument &e) {
        std::cerr << "minsup is not a positive integer: " << e.what()
            << std::endl;
        return EXIT_FAILURE;
    } catch (const std::out_of_range &e) {
        std::cerr << "minsup is not a positive integer: " << e.what()
            << std::endl;
        return EXIT_FAILURE;
    }
    if (fpgrowth::THRESHOLD <= 0) {
        std::cerr << "minsup is not a positive integer" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream ds {dsf};
    if (!ds.is_open()) {
        std::cerr << dsf << " could not be opened!" << std::endl;
        return EXIT_FAILURE;
	}
	int i;
	fpgrowth::FI_tree* fptree;


	fpgrowth::fp_buf=new fpgrowth::memory(60, 4194304L, 8388608L, 2);
	fptree = (fpgrowth::FI_tree*)fpgrowth::fp_buf->newbuf(1,
			sizeof(fpgrowth::FI_tree));
	fptree->init(-1, 0);
	fptree->scan1_DB(ds);
    ds.clear();
    ds.seekg(0, ds.beg);
	fpgrowth::ITlen = new int[fptree->itemno];
	fpgrowth::bran = new int[fptree->itemno];
	fpgrowth::compact = new int[fptree->itemno];
	fpgrowth::prefix = new int[fptree->itemno];

#ifdef CFI
	fpgrowth::list=new fpgrowth::stack(fptree->itemno, true);
#else
	fpgrowth::list=new fpgrowth::stack(fptree->itemno);
#endif

	assert(fpgrowth::list!=NULL && fpgrowth::bran!=NULL && fpgrowth::compact!=NULL && fpgrowth::ITlen!=NULL && fpgrowth::prefix!=NULL);

	for(i =0; i < fptree->itemno; i++)
	{
		fpgrowth::ITlen[i] = 0L;
		fpgrowth::bran[i] = 0;
	}

	fptree->scan2_DB(ds);


    //print the count of emptyset
#ifdef FI
    ftor(0, NULL, fpgrowth::TRANSACTION_NO);
#endif

#ifdef CFI
    if(fpgrowth::TRANSACTION_NO != fptree->count[0])
        ftor(0, NULL, fpgrowth::TRANSACTION_NO);
#endif


	if(fptree->Single_path())
	{
		fpgrowth::Fnode* node;
		int i=0;
		for(node=fptree->Root->leftchild; node!=NULL; node=node->leftchild)
		{
			fpgrowth::list->FS[i++]=node->itemname;
#ifdef CFI
			fpgrowth::list->counts[i-1] = node->count;
#endif
		}

#ifdef FI
        fptree->generate_all(fptree->itemno, ftor);
#endif

#ifdef CFI
        int Count;
        i=0;
        while(i<fptree->itemno)
        {
            Count = fpgrowth::list->counts[i];
            for(; i<fptree->itemno && fpgrowth::list->counts[i]==Count; i++);
            fpgrowth::ITlen[i-1]++;
            ftor(i, fpgrowth::list->FS, Count);
        }
#endif

#ifdef MFI
        ftor(fptree->itemno, fpgrowth::list->FS, fptree->head[fptree->itemno-1]->count);
        fpgrowth::ITlen[i-1]=1;
#endif
        delete fpgrowth::list;
		return EXIT_SUCCESS;
	}

	fpgrowth::current_fi = new bool[fptree->itemno];
	fpgrowth::supp=new int[fptree->itemno];		//for keeping support of items
	assert(fpgrowth::supp!=NULL&&fpgrowth::current_fi!=NULL);

	for(i = 0; i<fptree->itemno; i++)
	{
		fpgrowth::current_fi[i] = false;
		fpgrowth::supp[i]=0;
	}

#ifdef MFI
	fpgrowth::MFI_tree* LMFI;
	fpgrowth::mfitrees = (fpgrowth::MFI_tree**)new fpgrowth::MFI_tree*[fptree->itemno];
	fpgrowth::memory* Max_buf=new fpgrowth::memory(40, 1048576L, 5242880, 2);
		LMFI = (fpgrowth::MFI_tree*)Max_buf->newbuf(1, sizeof(fpgrowth::MFI_tree));
		LMFI->init(Max_buf, fptree, NULL, NULL, -1);
		fptree->set_max_tree(LMFI);
		fpgrowth::mfitrees[0] = LMFI;
		fptree->FPmax(ftor);
#endif

#ifdef CFI
		fpgrowth::CFI_tree* LClose;
		fpgrowth::cfitrees = (fpgrowth::CFI_tree**)new fpgrowth::CFI_tree*[fptree->itemno];
		fpgrowth::memory* Close_buf=new fpgrowth::memory(40, 1048576L, 5242880, 2);
		LClose = (fpgrowth::CFI_tree*)Close_buf->newbuf(1, sizeof(fpgrowth::CFI_tree));
		LClose->init(Close_buf, fptree, NULL, NULL, -1);
		fptree->set_close_tree(LClose);
		fpgrowth::cfitrees[0] = LClose;
		fptree->FPclose(ftor);
#endif

#ifdef FI
		fptree->FP_growth(ftor);
#endif

//	delete fp_buf;
	delete fpgrowth::list;
	delete []fpgrowth::current_fi;
	delete []fpgrowth::supp;

	return EXIT_SUCCESS;
}

#endif // _MAIN_H_
