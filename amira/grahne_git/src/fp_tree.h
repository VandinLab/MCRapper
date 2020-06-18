/*
   Author:  Jianfei Zhu
            Concordia University
   Date:    Sep. 26, 2003

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

#ifndef _FP_TREE_CLASS
#define _FP_TREE_CLASS

#include <fstream>
#include <vector>

#include "buffer.h"
#include "common.h"
#include "fp_node.h"

#define SORTHRESH 9

namespace fpgrowth {

class FP_tree {
public:
	int itemno;		//Header_table
	int *order;		//Header_table
	int *table;		//orders[table[i]]=i; table[orders[i]]=i;

public:
	~FP_tree(){}
};

class FI_tree : public FP_tree{
public:
	Fnode* Root;
	Fnode** head;	//Header_table

	int** array;   //the length of array is itemno-1; for each array[i], array[i][j] keeps the support of ij.
	int* count;		//Header_table
	double level_25;
	MFI_tree* LMaxsets;
	CFI_tree* LClose;

private:
	void insert(int* compact, int counts, int current);
	void scan1_DB(FI_tree*);			//build header_table
	void scan2_DB(FI_tree*, Fnode*);	//construct fp-tree
	int conditional_pattern_base(int, bool=false)const;
	int conditional_pattern_base(Fnode*, bool=false)const;
	memory* allocate_buf(int sequence, int iteration, int ms);
	void fill_count(int*, int);
	void cal_level_25();
	template<typename Ftor> void powerset(int*, int, int*, int, int, Ftor &)const;

public:
	void init(int Itemno, int new_item_no);
	void set_max_tree(MFI_tree* lmfisets){LMaxsets = lmfisets;}
	void set_close_tree(CFI_tree* lclose){LClose = lclose;}
	~FI_tree(){/*delete root;	delete []order;	delete []table;*/}

	void scan1_DB(std::ifstream &);		//find the count of all items from the dataset, and initialize a few data structures using scan1_DB_inity_arrays
    void scan1_DB_init_arrays(const std::vector<int> &, const int); // initialize a few data structures
	void scan2_DB(std:: ifstream &ds);		//construct the first fp-tree from  original DB
    //void FI_tree::scan2_DB(const std::vector<std::vector<int>>&);
    template<typename V> void scan2_DB(const V &ts);
	bool Single_path(bool=false)const;   //Is it a single path?
	template<typename Ftor> void generate_all(int, Ftor &)const;

	template<typename Ftor> int FP_growth(Ftor &ftor);
	template<typename Ftor> int FPmax(Ftor &ftor);
	template<typename Ftor> int FPclose(Ftor &ftor);
	//int FP_growth(FSout*);
	//int FPmax(FSout*);
	//int FPclose(FSout*);
};

class MFI_tree : public FP_tree{
public:
	Mnode* Root;
	Mnode **head;

	int MFSNo;
	int posi;                    //The position in list
	memory* Max_buf;
public:
	void init(memory*, FI_tree*, MFI_tree*, Mnode*, int);
	~MFI_tree();

	bool is_subset();
	void insert(bool* origin, int);
	void insert(int*, int, int);  //compact info
};

class CFI_tree : public FP_tree {
public:
	Cnode* Root;
	Cnode **head;

	memory* Close_buf;
	int CloseNo;
	int posi;                    //The position in list
public:
	void init(memory*, FI_tree*, CFI_tree*, Cnode*, int);
	~CFI_tree();

	bool is_subset(int Count)const;
	void insert(bool*, int, int);
	void insert(int*, int, int, int);  //compact info
	template<typename Ftor> bool generate_close(int, int, Ftor &);
	void order_FS(int*, int, int);
};

class stack{
public:
	int top;
	int* FS;
	int* counts;
public:
	stack(int, bool=false);
	~stack();

	void insert(FI_tree* fptree);
};

} // namespace fpgrowth

#endif
