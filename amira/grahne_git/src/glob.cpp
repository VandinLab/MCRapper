
#include "common.h"

namespace fpgrowth {

int *ITlen;
int* bran;
int* prefix;

int* order_item;		// given order i, order_item[i] gives itemname
int* item_order;		// given item i, item_order[i] gives its new order
						//	order_item[item_order[i]]=i; item_order[order_item[i]]=i;
bool* current_fi;
int* compact;
int* supp;

MFI_tree** mfitrees;
CFI_tree** cfitrees;

stack* list;
int TRANSACTION_NO=0;
int ITEM_NO=100;
int THRESHOLD;

memory* fp_buf;

} // namespace fpgrowth
