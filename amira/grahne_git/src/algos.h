#ifndef _ALGOS_H_
#define _ALGOS_H_

#include <vector>

#include "common.h"
#include "fp_tree.h"

namespace fpgrowth {

template<typename T> void FI_tree::scan2_DB(const T &ts) {
    int* origin;
    origin = new int[ITEM_NO];
    for(int j=0; j<ITEM_NO; j++)
        origin[j]=-1;


    for (const auto &p : ts) {
        if constexpr (std::is_same_v<T,
                std::vector<std::vector<unsigned int>>>) {
            // T is a vector of transactions.
            int has {0};
            for (const auto i : p)  {
                if (item_order[i] != -1) {
                    ++has;
                    origin[item_order[i]] = item_order[i];
                }
            }
            if(has){
                fill_count(origin, 1);
                insert(compact, 1, has);
            }
        } else {
            // T is a map of unique transactions and their counts
            for (unsigned int j {0}; j < p.second; ++j) {
                int has {0};
                for (const auto i : p.first)  {
                    if (item_order[i] != -1) {
                        ++has;
                        origin[item_order[i]] = item_order[i];
                    }
                }
                if(has){
                    fill_count(origin, 1);
                    insert(compact, 1, has);
                }
            }
        }
    }

    cal_level_25();

    delete []origin;
}

void update_mfi_trees(int treeno)
{
	int i, j;
	for(i=0; i<treeno; i++)
	{
		for(j = list->top-1; j >mfitrees[i]->posi; j--)
			current_fi[mfitrees[i]->order[list->FS[j]]] = true;

		mfitrees[i]->insert(current_fi, list->top-1-mfitrees[i]->posi);
	}
}

void update_cfi_trees(int treeno, int Count)
{
	int i, j;
	for(i=0; i<treeno; i++)
	{
		for(j = list->top-1; j >cfitrees[i]->posi; j--)
			current_fi[cfitrees[i]->order[list->FS[j]]] = true;

		cfitrees[i]->insert(current_fi, Count, list->top-1-cfitrees[i]->posi);
	}
}

template<typename Ftor> void FI_tree::powerset(int* prefix, int prefixlen,
        int* items, int current, int itlen, Ftor &ftor) const {
    if (current == itlen) {
        if(prefixlen != 0) {
            ITlen[list->top+prefixlen-1]++;
            ftor(prefixlen, prefix,
                    this->head[order[prefix[prefixlen-1]]]->count, list->top,
                    list->FS);
        }
    } else {
        current++;
        powerset(prefix, prefixlen, items, current, itlen, ftor);
        current--;
        prefix[prefixlen++] = items[current++];
        powerset(prefix, prefixlen, items, current, itlen, ftor);
    }
}

template<typename Ftor> void FI_tree::generate_all(int new_item_no, Ftor &ftor) const {
	powerset(prefix, 0, list->FS, list->top, list->top+new_item_no, ftor);
}

template<typename Ftor> bool CFI_tree::generate_close(int new_item_no,
        int Count, Ftor &ftor) {
	int i, whole, temp = Count;

	if(list->top>1)
		temp=list->counts[list->top-2];

	whole = list->top+new_item_no;
	for(i=list->top; i<whole && list->counts[i]==Count; i++);
	list->top = i;

	ITlen[i-1]++;
    ftor(i, list->FS, Count);

//	insert(list->FS, posi+1, i, Count);
//	update_cfi_trees(posi+1, Count);
	update_cfi_trees(posi+2, Count);  //Modified on Oct. 10, 2003

	if(i==whole && temp == Count)
		return true;

	while(i<whole)
	{
		Count = list->counts[i];
		for(; i<whole && list->counts[i]==Count; i++);
		list->top = i;

//Patch Dec. 4
		order_FS(list->FS, posi+2, i-1);
		if(!is_subset(Count))
		{
			ITlen[i-1]++;
			ftor(i, list->FS, Count);
//			insert(list->FS, posi+1, i, Count);
//			update_cfi_trees(posi+1, Count);
			update_cfi_trees(posi+2, Count);   //Modified on Oct. 10, 2003
		}else{
		}
	}

	return false;
}


//int FI_tree::FPmax(FSout* fout)
template<typename Ftor> int FI_tree::FPmax(Ftor &ftor)
{
	static int ms=9;		//power2[i] is the smallest block size for memory  2**9 = 512

	int i, sequence, new_item_no, listlen;
	int MC=0;			//markcount for memory
	unsigned int MR=0;	//markrest for memory
	char* MB;			//markbuf for memory

	Fnode* Current;

	for(sequence=itemno-1; sequence>=0; sequence--)
	{
		Current=head[sequence];
		list->FS[list->top++]=head[sequence]->itemname;
		listlen = list->top;

		if(array && sequence>SUDDEN+1)
			new_item_no=conditional_pattern_base(Current->itemname);  //new_item_no is the number of elements in new header table.
		else if(sequence !=0)
			new_item_no=conditional_pattern_base(Current);  //new_item_no is the number of elements in new header table.
		else
			new_item_no = 0;


		if(LMaxsets->is_subset())
		{
			for(int j=listlen; j < list->top; j++)supp[order[list->FS[j]]] = 0;

			list->top=listlen-1;
			if(new_item_no == sequence)return new_item_no;
			continue;
		}
		if(new_item_no==0 || new_item_no == 1)
		{
			ITlen[list->top-1]++;
            if(new_item_no==1)
                ftor(list->top, list->FS, supp[order[list->FS[list->top-1]]]);
            else
                //fout->printSet(list->top, list->FS, this->head[sequence]->count);
                ftor(list->top, list->FS, this->head[sequence]->count);
			LMaxsets->insert(list->FS, LMaxsets->posi+1, list->top+new_item_no-1);
			update_mfi_trees(LMaxsets->posi+1);
			list->top = listlen+1;
			if(new_item_no==1)supp[order[list->FS[list->top-1]]] = 0;
			list->top=listlen-1;
			if(new_item_no != sequence)continue;
			return new_item_no;
		}

		FI_tree *fptree;
		MB=fp_buf->bufmark(&MR, &MC);

		fptree = (FI_tree*)fp_buf->newbuf(1, sizeof(FI_tree));
		fptree->init(this->itemno, new_item_no);

		fptree->scan1_DB(this);
		fptree->scan2_DB(this, Current);

		list->top=listlen;
		if(fptree->Single_path())
		{
            ftor(list->top+new_item_no, list->FS, fptree->head[fptree->itemno-1]->count);
			ITlen[list->top+new_item_no-1]++;
			LMaxsets->insert(list->FS, LMaxsets->posi+1, list->top+new_item_no);
			list->top = list->top+new_item_no;
			update_mfi_trees(LMaxsets->posi+1);
			list->top = listlen;
			list->top--;
			if(new_item_no == sequence)
			{
				fp_buf->freebuf(MR, MC, MB);
				return new_item_no;
			}
		}else{             //not single path
			memory* Max_buf;
			Max_buf=allocate_buf(sequence, LMaxsets->posi, ms);

			MFI_tree* new_LMFI = (MFI_tree*)Max_buf->newbuf(1, sizeof(MFI_tree));
			new_LMFI->init(Max_buf, fptree, LMaxsets, LMaxsets->head[sequence], list->top-1);
			fptree->set_max_tree(new_LMFI);
			mfitrees[LMaxsets->posi+2] = new_LMFI;
			i=fptree->FPmax(ftor);

			list->top = new_LMFI->posi;
			if(Max_buf->half())ms++;
			delete Max_buf;

			if(i+1 == sequence)
			{
				fp_buf->freebuf(MR, MC, MB);
				return i+1;
			}
		}
		fp_buf->freebuf(MR, MC, MB);
	}
	return 0;
}

//int FI_tree::FPclose(FSout* fout)
template<typename Ftor> int FI_tree::FPclose(Ftor &ftor) {
	static int ms=9;		//power2[i] is the smallest block size for memory  2**9 = 512

	int sequence, new_item_no, listlen;
	int MC=0;			//markcount for memory
	unsigned int MR=0;	//markrest for memory
	char* MB;			//markbuf for memory

	Fnode* Current;

	for(sequence=itemno-1; sequence>=0; sequence--)
	{
		Current=head[sequence];
		list->FS[list->top++]=head[sequence]->itemname;
		list->counts[list->top-1] = this->head[sequence]->count;
		listlen = list->top;

		if(LClose->is_subset(this->head[sequence]->count))
		{
			list->top=listlen-1;
			continue;
		}
		if(array && sequence>SUDDEN+1)
			new_item_no=conditional_pattern_base(Current->itemname, 1);  //new_item_no is the number of elements in new header table.
		else if(sequence !=0)
			new_item_no=conditional_pattern_base(Current, 1);  //new_item_no is the number of elements in new header table.
		else
			new_item_no = 0;

		if(new_item_no==0 || new_item_no == 1)
		{
			if(new_item_no==0)
			{
				ITlen[list->top-1]++;
                ftor(list->top, list->FS, this->head[sequence]->count);
				LClose->insert(list->FS, LClose->posi+1, list->top+new_item_no-1, this->head[sequence]->count);
				update_cfi_trees(LClose->posi+1, this->head[sequence]->count);
			}else{
				list->top=listlen;
//Bug! Dec. 5
/*				if(LClose->generate_close(1, this->head[sequence]->count, fout))
				{
					supp[order[list->FS[listlen]]] = 0;
					list->top=listlen-1;
					return new_item_no;
				}else
					supp[order[list->FS[listlen]]] = 0;
*/				LClose->generate_close(1, this->head[sequence]->count, ftor);
				supp[order[list->FS[listlen]]] = 0;
			}
			list->top=listlen-1;
			continue;
		}
/*		int j;
		for(j=listlen; j<new_item_no+listlen && list->counts[listlen-1]!=list->counts[j]; j++);
		if(j==new_item_no+listlen)
		{
			ITlen[listlen-1]++;
			if(fout)
				fout->printSet(listlen, list->FS, this->head[sequence]->count);

			LClose->insert(list->FS, LClose->posi+1, listlen, this->head[sequence]->count);
			list->top = listlen;
			update_cfi_trees(LClose->posi+1, this->head[sequence]->count);
			list->top += new_item_no;
		}
*/
		FI_tree *fptree;
		MB=fp_buf->bufmark(&MR, &MC);

		fptree = (FI_tree*)fp_buf->newbuf(1, sizeof(FI_tree));
		fptree->init(this->itemno, new_item_no);

		fptree->scan1_DB(this);
		fptree->scan2_DB(this, Current);

		list->top=listlen;
		if(fptree->Single_path(1))
		{
			list->top=listlen;
			if(LClose->generate_close(new_item_no, this->head[sequence]->count, ftor))
			{
				if(new_item_no == sequence)
				{
					fp_buf->freebuf(MR, MC, MB);
					list->top = listlen-1;
					return new_item_no;
				}
			}
			list->top = listlen-1;
		}else{             //not single path
                               /* patch   Oct. 9, 2003 */
                        int j;
                        for(j=listlen; j<new_item_no+listlen && list->counts[listlen-1]!=list->counts[j]; j++);
                        if(j==new_item_no+listlen)
                        {
                                ITlen[listlen-1]++;
                                ftor(listlen, list->FS, this->head[sequence]->count);
                                LClose->insert(list->FS, LClose->posi+1, listlen, this->head[sequence]->count);
                                list->top = listlen;
                                update_cfi_trees(LClose->posi+1, this->head[sequence]->count);
//                              list->top += new_item_no;
                        }
                              /***************************/
			memory* Close_buf;
			Close_buf=allocate_buf(sequence, LClose->posi, ms);

			CFI_tree* new_LClose = (CFI_tree*)Close_buf->newbuf(1, sizeof(CFI_tree));
			new_LClose->init(Close_buf, fptree, LClose, LClose->head[sequence], list->top-1);
			fptree->set_close_tree(new_LClose);
			cfitrees[LClose->posi+2] = new_LClose;
			fptree->FPclose(ftor);

			list->top = new_LClose->posi;
			if(Close_buf->half())ms++;
			delete Close_buf;

		}
		fp_buf->freebuf(MR, MC, MB);
	}
	return 0;
}

//int FI_tree::FP_growth(FSout* fout)
template<typename Ftor> int FI_tree::FP_growth(Ftor &ftor) {
	int sequence, new_item_no, listlen;
	int MC=0;			//markcount for memory
	unsigned int MR=0;	//markrest for memory
	char* MB;			//markbuf for memory

	Fnode* Current;

	for(sequence=itemno-1; sequence>=0; sequence--)
	{
		Current=head[sequence];
		list->FS[list->top++]=head[sequence]->itemname;
		listlen = list->top;

		ITlen[list->top-1]++;
        ftor(list->top, list->FS, this->head[sequence]->count);

		if(array && sequence>SUDDEN+1)
			new_item_no=conditional_pattern_base(Current->itemname);  //new_item_no is the number of elements in new header table.
		else if(sequence !=0)
			new_item_no=conditional_pattern_base(Current);  //new_item_no is the number of elements in new header table.
		else
			new_item_no = 0;

		if(new_item_no==0 || new_item_no == 1)
		{
			if(new_item_no==1)
			{
				ITlen[list->top-1]++;
                ftor(list->top, list->FS, supp[order[list->FS[list->top-1]]]);
			}
			if(new_item_no==1)supp[order[list->FS[list->top-1]]] = 0;
			list->top=listlen-1;
			continue;
		}

		FI_tree *fptree;
		MB=fp_buf->bufmark(&MR, &MC);

		fptree = (FI_tree*)fp_buf->newbuf(1, sizeof(FI_tree));
		fptree->init(this->itemno, new_item_no);

		fptree->scan1_DB(this);
		fptree->scan2_DB(this, Current);

		list->top=listlen;
		if(fptree->Single_path())
		{
                               /* patch   Oct. 9, 2003*/
                        Fnode* node;
                        for(node=fptree->Root->leftchild; node!=NULL; node=node->leftchild)
                                list->FS[list->top++] = node->itemname;
                        list->top = listlen;
                              /*************************/
			fptree->generate_all(new_item_no, ftor);
			list->top--;
		}else{             //not single path
			fptree->FP_growth(ftor);
			list->top = listlen-1;
		}
		fp_buf->freebuf(MR, MC, MB);
	}
	return 0;
}


} // namespace fpgrowth

#endif // _ALGOS_H_
