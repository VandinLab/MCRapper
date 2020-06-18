/* Linear time Closed itemset Miner for Frequent Itemset Mining problems */
/* 2004/4/10 Takeaki Uno */
/* This program is available for only academic use.
   Neither commercial use, modification, nor re-distribution is allowed */

#ifndef _lcm_c_
#define _lcm_c_

#include<time.h>
#include"lib_e.c"
#include"lcm_var.c"
#define LCM_PROBLEM LCM_CLOSED
#include"trsact.c"
#include"lcm_io.c"
#include"lcm_init.c"
#include"lcm_lib.c"
#include <unistd.h>

/* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
#include"transaction_keeping.c"
/* END OF MODIFICATIONS */

/* MODIFICATIONS FOR FAST WY ALGORITHIM */
//#include"wy.c"
#include"mcrapper.c"
#include"permutation.c"
/* END OF MODIFICATIONS */

/* MODIFICATIONS TO KEEP TRACK OF EXECUTION TIME AND MEMORY CONSUMPTION */
#include"time_keeping.c"
/* END OF MODIFICATIONS */

/* FUNCTION DECLARATIONS OF ORIGINAL LCM SOURCE */
void LCMclosed_BM_iter(int item, int m, int pmask, struct bounds_info* parent_info , struct bounds_info* current_info);

void LCMclosed_BM_recursive (int item, int mask, int pmask  , struct bounds_info* parent_info , struct bounds_info* current_info){
  /*printf("Call to LCMclosed_BM_recursive. item %d\n",item);

  printf("  parent info supp %d\n",parent_info->frequency);
  printf("  parent info as %d\n",parent_info->as);
  printf("  parent info max_as %d\n",parent_info->max_as);
  printf("  parent info min_as %d\n\n",parent_info->min_as);
  printf("  current info supp %d\n",current_info->frequency);
  printf("  current info as %d\n",current_info->as);
  printf("  current info max_as %d\n",current_info->max_as);
  printf("  current info min_as %d\n\n",current_info->min_as);*/
  int i;
  if((i=LCM_Ofrq[0]) >= LCM_th){
    if((LCM_BM_pp[1]&pmask)==1 && LCM_BM_pt[1]==0){
        struct bounds_info child_info;
        init_bounds(&child_info);
        child_info.frequency = LCM_Ofrq[0];
        reset_bounds(&child_info);
        //printf("calling LCM_print_last with supp %d\n", LCM_Ofrq[0]);
        LCM_print_last(LCM_Op[0], i , current_info , &child_info);
    }
//     else printf ("11clo item%d %d: %d : %x & %x = %x,   %d:  prv%d pprv%d \n", LCM_Op[0], LCM_itemsett, LCM_frq, pmask, LCM_BM_pp[1], LCM_BM_pp[1] & pmask, LCM_BM_pt[1], LCM_prv, LCM_pprv);
   // pruning has to be here
  }
  LCM_BM_weight[1] = LCM_Ofrq[0] = 0;
  LCM_Ot[0] = LCM_Os[0];
  /* MODIFICATIONS TO KEEP TRACK OF TRANSACTIONS */
  BM_TRANS_LIST_EMPTY(1);
  /* END OF MODIFICATIONS */
  for(i=1; i<item; i++)
    if(LCM_Ofrq[i] >= LCM_th) {
      struct bounds_info child_info;
      init_bounds(&child_info);
      child_info.frequency = LCM_Ofrq[i];
      reset_bounds(&child_info);
      //printf("   calling LCMclosed_BM_iter from LCMclosed_BM_recursive with item %d\n", i);
      LCMclosed_BM_iter(i, mask, pmask , current_info , &child_info);
      child_info.frequency = -1; // debug
    }
    else{
    	if(LCM_Ofrq[i]>0) {
    		//printf("Item %d is infrequent (%d,%d) and will not have LCMclosed_BM_iter,\n",i,LCM_Ofrq[i],LCM_th);
    		LCM_BM_occurrence_delete(i);
    		//printf("Problem fixed!\n");
    	}
    }
  	 // Maybe an else if (LCM_Ofrq[i] > 0) LCM_BM_occurrence_delete(i) doesn't hurt
     //printf("Call to LCMclosed_BM_recursive done. item %d\n",i);
}



/*************************************************************************/
/* LCMclosed iteration (bitmap version ) */
/* input: T:transactions(database), item:tail(current solution) */
/*************************************************************************/
void LCMclosed_BM_iter(int item, int m, int pmask , struct bounds_info* parent_info , struct bounds_info* current_info){
  /*printf("Call to LCMclosed_BM_iter. item %d\n",item);

  printf("  parent info supp %d\n",parent_info->frequency);
  printf("  parent info as %d\n",parent_info->as);
  printf("  parent info max_as %d\n",parent_info->max_as);
  printf("  parent info min_as %d\n\n",parent_info->min_as);
  printf("  current info supp %d\n",current_info->frequency);
  printf("  current info as %d\n",current_info->as);
  printf("  current info max_as %d\n",current_info->max_as);
  printf("  current info min_as %d\n\n",current_info->min_as);*/

  int mask, it = LCM_itemsett, ttt;

  LCM_frq = LCM_Ofrq[item];

  //printf("  supp %d\n",LCM_frq);

  if(LCM_frq != current_info->frequency){
    printf("problem in LCMclosed_BM_iter: supp %d, current supp %d\n",LCM_frq,current_info->frequency);
    printf("  parent info supp %d\n",parent_info->frequency);
    //printf("  parent info as %d\n",parent_info->as);
    printf("  parent info max_as %d\n",parent_info->max_as);
    //printf("  parent info min_as %d\n\n",parent_info->min_as);
    printf("  current info supp %d\n",current_info->frequency);
    //printf("  current info as %d\n",current_info->as);
    printf("  current info max_as %d\n",current_info->max_as);
    //printf("  current info min_as %d\n\n",current_info->min_as);
    current_info->frequency = LCM_frq;
    reset_bounds(current_info);
  }

  pmask &= BITMASK_31[item];
  if((ttt = LCM_BM_closure(item, pmask)) > 0){
   // pruning has to be here
	  //printf ("BMclo %d item%d it%d frq%d,  prv%d pprv%d::  ttt=%d,%d pmask%x\n", item, LCM_Op[item], LCM_itemsett, LCM_frq, LCM_prv, LCM_pprv, ttt,LCM_Op[ttt], pmask );
    LCM_BM_occurrence_delete(item);
    return;
  }
  LCM_iters++;
  BUF_reset(&LCM_B);
  LCMclosed_BM_occurrence_deliver_(item, m);
  LCM_additem(LCM_Op[item]);
  mask = LCM_BM_rm_infreq(item, &pmask);

  LCM_solution();

  /* MODIFICATION FOR FAST WY ALGORITHIM */
  if(LCM_frq != current_trans.siz){
	  printf("LCM_frq=%d, current_trans.siz=%d\n",LCM_frq,current_trans.siz);
  }
  bm_process_solution(LCM_frq,item,&mask,parent_info,current_info);
  /* END OF MODIFICATION */

  /* MODIFICATION TO KEEP TRACK OF TRANSACTIONS */
  //print_current_trans();
  BM_CURRENT_TRANS_EMPTY();
  /* END OF MODIFICATION */

  LCMclosed_BM_recursive(item, mask, pmask,parent_info,current_info);
  while(LCM_itemsett>it) LCM_delitem();
  BUF_clear(&LCM_B);

  //printf("Call to LCMclosed_BM_iter done. item %d\n",item);
}




/***************************************************************/
/* iteration of LCM ver. 2 */
/* INPUT: T:transactions(database), item:tail of the current solution */
/*************************************************************************/
// LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
int LCMclosed_iter(ARY *T, int item, int prv, TRANS_LIST *trans_list , struct bounds_info* parent_info , struct bounds_info* current_info){
  ARY TT;
  int i, ii, e, ee, n, js=LCM_jump.s, it=LCM_itemsett, mask;
  int flag=-1, perm[LCM_BM_MAXITEM], pmask = 0xffffffff;
  QUEUE_INT *q;

  /*printf("Call to LCMclosed_iter. item %d\n",item);
  printf("  parent info supp %d\n",parent_info->frequency);
  printf("  parent info as %d\n",parent_info->as);
  printf("  parent info max_as %d\n",parent_info->max_as);
  printf("  parent info min_as %d\n\n",parent_info->min_as);
  printf("  current info supp %d\n",current_info->frequency);
  printf("  current info as %d\n",current_info->as);
  printf("  current info max_as %d\n",current_info->max_as);
  printf("  current info min_as %d\n\n",current_info->min_as);*/

  LCM_jump.s = LCM_jump.t;
  LCM_iters++;
  LCM_additem(item);
  LCM_frq = LCM_Ofrq_[item];
  LCM_prv = item;
  LCM_pprv = prv;

  //printf("   supp %d\n",LCM_frq);

  if(LCM_frq != current_info->frequency && current_info->frequency > 0){
    printf("problem in LCMclosed_iter: supp %d, current supp %d\n",LCM_frq,current_info->frequency);
    printf("Call to LCMclosed_iter. item %d\n",item);
    printf("  parent info supp %d\n",parent_info->frequency);
    //printf("  parent info as %d\n",parent_info->as);
    printf("  parent info max_as %d\n",parent_info->max_as);
    //printf("  parent info min_as %d\n\n",parent_info->min_as);
    printf("  current info supp %d\n",current_info->frequency);
    //printf("  current info as %d\n",current_info->as);
    printf("  current info max_as %d\n",current_info->max_as);
    //printf("  current info min_as %d\n\n",current_info->min_as);
    current_info->frequency = LCM_frq;
    reset_bounds(current_info);
  }

  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  TRANS_LIST mk_trans_list, shrink_trans_list;
  /* END OF MODIFICATIONS */

  //printf( " Ot %d %d (%d,%d)\n", LCM_Ot[item]-LCM_Os[item], LCM_frq, item,prv);
  n = LCM_freq_calc(T, item, LCM_Eend-1);
  if(prv >= 0) LCM_Ofrq[prv] = 0;
  LCM_Ofrq[item] = 0;
  ii = LCM_jump_rm_infreq(item);
  LCM_jumpt = LCM_jump.t;
  if(ii > item){
    flag = ii;
    //printf ("###clo item%d %d %d: %d\n", item, ii, LCM_Ofrq[ii], LCM_Ofrq_[ii]);
    goto END2;
  }  /* itemset is not closed */

  BUF_reset(&LCM_B);
  LCM_partition_prefix(item);

  if(QUEUE_LENGTH(LCM_jump)==0){
    LCM_Ofrq[item] = 0;
    LCMclosed_BM_occurrence_deliver_first(item, T, trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
    for(i=LCM_jump.s; i<LCM_jumpt; i++) LCM_Ofrq[LCM_jump.q[i]] = 0;
    mask = LCM_BM_rm_infreq(LCM_BM_MAXITEM, &pmask);
    LCM_solution();
    /* MODIFICATIONS FOR WY ALGORITHM */
    ary_process_solution(LCM_frq, trans_list, item, &mask, parent_info , current_info);
    /* END OF MODIFICATIONS */
    /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
    //print_transaction_list(trans_list,item);
    /* END OF MODIFICATIONS */
    struct bounds_info child_info__;
    init_bounds(&child_info__);
    LCMclosed_BM_recursive(LCM_BM_MAXITEM, mask, pmask, current_info, &child_info__ );
    BUF_clear(&LCM_B);
    goto END2;
  }

  LCM_BM_occurrence_deliver_first(item, T);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
  QUEUE_FE_LOOP_(LCM_jump, i, ii) LCM_Ofrq[ii] = 0;
  mask = LCM_BM_rm_infreq(LCM_BM_MAXITEM, &pmask);
  LCM_solution();
  /* MODIFICATIONS FOR WY ALGORITHM */
  ary_process_solution(LCM_frq, trans_list, item, &mask, parent_info , current_info);
  /* END OF MODIFICATIONS */
  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  //print_transaction_list(trans_list,item);
  /* END OF MODIFICATIONS */
  BUF_clear(&LCM_B);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  TRANS_LIST_INIT(&mk_trans_list, LCM_frq, LCM_Ot[item]-LCM_Os[item]);
  /* END OF MODIFICATIONS */
  for(i=0; i<LCM_BM_MAXITEM; i++) perm[i] = LCM_Op[i];
  QUEUE_FE_LOOP_(LCM_jump, i, ii) LCM_Ofrq[ii] = LCM_th;
  LCM_Ofrq[item] = LCM_th;
  /* LAST TWO ARGUMENTS ADDED FOR TRANSACTION KEEPING */
  LCM_mk_freq_trsact(&TT, T, item, LCM_Eend-1, n+(LCM_Ot[item]-LCM_Os[item]), mask, trans_list, &mk_trans_list);
  /* END OF MODIFICATIONS */
  LCM_Ofrq[item] = 0;

  BUF_reset(&LCM_B);
  for(i=0; i<LCM_BM_MAXITEM; i++) LCM_BM_occurrence_delete(i);
  LCMclosed_BM_occurrence_deliver_first(-1, &TT, &mk_trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
  for(i=LCM_jump.s; i<LCM_jumpt; i++) LCM_Ofrq[LCM_jump.q[i]] = 0;
  struct bounds_info child_info_;
  init_bounds(&child_info_);
  LCMclosed_BM_recursive(LCM_BM_MAXITEM, 0xffffffff, 0xffffffff, current_info, &child_info_ );
  BUF_clear(&LCM_B);

  if(QUEUE_LENGTH(LCM_jump) == 0) goto END0;
  q = ((QUEUE *)(TT.h))->q;
  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  if(ii >= 2 && TT.num>5){
	  TRANS_LIST_INIT(&shrink_trans_list, mk_trans_list.siz1, mk_trans_list.siz2);
	  LCM_shrink(&TT, item, 1, &mk_trans_list, &shrink_trans_list);//LAST TWO ARGUMENTS ADDED FOR TRANSACTION KEEPING
	  TRANS_LIST_END(&mk_trans_list);
  }else{
	  shrink_trans_list = mk_trans_list;
  }
  /* END OF MODIFICATIONS */
  LCM_occ_deliver(&TT, item-1);

  do{
    i = QUEUE_ext_tail_(&LCM_jump);
//    printf ("i=%d(%d) %d %d  :%d,%d\n", i, item, LCM_Ot[item]-LCM_Os[item], LCM_Ot[i]-LCM_Os[i], LCM_Ofrq[i], LCM_Ofrq_[i]);
    /* MODIFICATIONS FOR WY ALGORITHM */
    //WY permutations might cause that some items in LCM_jump.q are no longer frequent/testable, so
    // a check must be added
    if(LCM_Ofrq_[i]>=LCM_th) {
      struct bounds_info child_info_i;
      init_bounds(&child_info_i);
      child_info_i.frequency = LCM_Ofrq_[i];
      reset_bounds(&child_info_i);
      //printf("calling LCMclosed_iter 1\n");
      //printf("  LCM_Ofrq_[i] %d\n",LCM_Ofrq_[i]);
      /*printf("  child_info_i.frequency %d\n",child_info_i.frequency);
      printf("  current_info->frequency %d\n\n",current_info->frequency);*/
      ii = LCMclosed_iter(&TT, i, item, &shrink_trans_list,current_info , &child_info_i);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
    }
    //else printf("Item %d pruned since LCM_Ofrq_[%d]=%d, and LCM_th=%d\n",i,i,LCM_Ofrq_[i],LCM_th);
    /* END OF MODIFICATIONS */
    LCM_Ot[i] = LCM_Os[i];
    LCM_Ofrq_[i] = 0;
  }while(LCM_jump.t > LCM_jump.s);

  free2(q);
  ARY_end(&TT);
  TRANS_LIST_END(&shrink_trans_list);
  END0:;
  for(i=0; i<LCM_BM_MAXITEM; i++) LCM_Op[i] = perm[i];
  goto END3;
  END2:;
  for(i=LCM_jump.s; i<LCM_jumpt; i++) LCM_Ofrq[LCM_jump.q[i]] = 0;
  LCM_jump.t = LCM_jump.s;
  END3:;
  LCM_jump.s = js;
  while(it<LCM_itemsett) LCM_delitem();

  //printf("Call to LCMclosed_iter done. item %d\n",item);


  return (flag);
}

/***************************************************************/
/* main of LCM ver. 3 */
/*************************************************************************/
void LCMclosed(){
  int i;
  BUF_reset(&LCM_B);
  LCMclosed_BM_occurrence_deliver_first(-1, &LCM_Trsact, &root_trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
  for (i=LCM_BM_MAXITEM ; i<LCM_Eend ; i++){
    LCM_Ofrq_[i] = LCM_Ofrq[i];
    LCM_Ofrq[i] = 0;
  }

  struct bounds_info child_info;
  init_bounds(&child_info);
  struct bounds_info current_info;
  init_bounds(&current_info);

  LCMclosed_BM_recursive(LCM_BM_MAXITEM, 0xffffffff, 0xffffffff, &current_info, &child_info);
  BUF_clear(&LCM_B);

  for (i=LCM_BM_MAXITEM ; i<LCM_Eend ; i++){
    struct bounds_info child_info_i;
    init_bounds(&child_info_i);
    //printf("calling LCMclosed_iter 2\n");
    LCMclosed_iter(&LCM_Trsact, i, -1, &root_trans_list, &current_info, &child_info_i);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
    LCM_Ot[i] = LCM_Os[i];
    LCM_Ofrq_[i] = LCM_Ofrq[i] = 0;
  }
  LCM_iters++;
}

int random_seed;
char *input_file; // dataset path
double max_eps;
char *basefilename;

int parseParameters(int argc, char *argv[]){

	// default values for parameters
	K = -1; // -1 means unbounded search

	delta = 0.1;
	J = 10; // number of draws of Rademacher random variables
  int input_file_given = 0;
  int base_file_given = 0;
  int print_help = 0;
  random_seed = 1;
  max_eps = 1.0;
  theta = 0.0;
  max_pattern_freq = 1.0;
  max_variance = 0.25;
  dev_lower_bound = 0.0;
  approximation_ratio = 1.0;
  union_bound_term = 0.0;

	char opt;
  while ( ( opt = getopt ( argc, argv, "j:d:b:i:k:s:e:h:t:f:v:l:a:u:" ) ) != -1 ) {
    switch ( opt ) {
    case 'j': J = atoi(optarg); break;
    case 'd': delta = atof(optarg); break;
    case 'i': input_file = optarg; input_file_given = 1; break;
    case 'b': basefilename = optarg; base_file_given = 1; break;
    case 'k': K = atoi(optarg); break;
		case 's': random_seed = atoi(optarg); break;
    case 'e': max_eps = atof(optarg); break;
    case 'h': print_help = 1; break;
    case 't': theta = atof(optarg); break;
    case 'f': max_pattern_freq = atof(optarg); break;
    case 'v': max_variance = atof(optarg); break;
    case 'l': dev_lower_bound = atof(optarg); break;
    case 'a': approximation_ratio = atof(optarg); break;
    case 'u': union_bound_term = atof(optarg); break;
    }
  }

	// check for the spec file, prompt the parameters info
	if((!input_file_given || !base_file_given) || print_help){
    printf("Parameters:\n");
    printf("    -i [input file]\n");
    printf("    -d [probability delta] (in (0,1] , default: 0.1)\n");
    printf("    -j [Number of Rademacher vectors] (default: 10)\n");
    printf("    -k [Number of top-k patterns to enumerate] (default: -1 (unbounded search))\n");
    printf("    -s [random seed] (default: 1)\n");
    printf("    -t [theta] (default: 1)\n");
    printf("    -e [maximum acceptable value of epsilon] (default: 1)\n");
    printf("    -f [maximum frequency of patterns to use for n-MCERA] (default: 1.)\n");
    printf("    -f [maximum variance] (default: 0.25)\n");
    printf("    -f [approximation ratio] (default: 1.0)\n");
    printf("    -l [lower bound to the value of the n-MCERA] (default: 0.0)\n");
    printf("    -h [print help]\n");
		return 1;
	}

	// check parameters correctness
	/*ifstream specfile (input_file);
	if (!specfile.is_open()){
		cout << "Unable to open specification file\n";
		return 0;
	}
	k=K_significant_patterns;
	k_max=k;
	K_significant_patterns_total = K_significant_patterns;
	if ( k < 1 && k != -1 ){
		cout << "Invalid value for k\n" << endl;
		return 0;
	}
	if ( jp < 1 ){
		cout << "Invalid value for jp!" << endl;
		return 0;
	}
	if ( alpha < 0.0 ){
		cout << "Invalid value for alpha\n" << endl;
		return 0;
	}
	if ( generalized_FWER < 1 ){
		cout << "Invalid value for generalized_FWER!" << endl;
		return 0;
	}
	if ( max_ram < 0 ){
			cout << "Invalid value for max ram\n" << endl;
			return 0;
	}*/

	/*cout << "dataset transactions = " << fileinput << endl;
	cout << "dataset labels = " << c_fileinput << endl;
	cout << "k = " << K_significant_patterns << endl;
	cout << "alpha = " << alpha << endl;*/
  printf("delta %f\n",delta);
  printf("number of rade draws %d\n",J);
  /*cout << "delta = " << delta << endl;
	cout << "jp = " << jp << endl;*/
	return 0;
}

/*************************************************************************/
/*************************************************************************/
int main(int argc, char *argv[]){
  int i;

  /* MODIFICATIONS FOR FAST WY ALGORITHIM */

  // Main input arguments which are not part of LCM
  //int n_perm;
  double sig_th;
  //char *class_labels_file;
  //char *trans_filename;
  char *tmp_filename;
  //int seed_idx;

  if(parseParameters(argc, argv) == 1){
    exit(1);
  }

  /*struct bounds_info test;
  printf("test: %d\n",test.frequency);*/


  // Initial time
  t_init = measureTime();

  // Check if input contains all needed arguments
  /*printf("argc %d\n",argc);
  if (argc < 5){
	  printf("LCM_WY_FISHER: output_basefilename n_draws delta transactions_file seed [k]\n");
	  exit(1);
  }*/

  // Create output files for results and profiling
  //tmp_filename = (char *)malloc((strlen(argv[1])+512)*sizeof(char));
  tmp_filename = (char *)malloc((strlen(basefilename)+512)*sizeof(char));
  if(!tmp_filename){
	fprintf(stderr,"Error in function main: couldn't allocate memory for array tmp_filename\n");
	exit(1);
  }
  // Create a file to report runtime information
  strcpy(tmp_filename,basefilename); strcat(tmp_filename,"_timing.txt");
  if(!(timing_file = fopen(tmp_filename,"w"))){
	fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
	exit(1);
  }
  // Create a file to report results
  strcpy(tmp_filename,basefilename); strcat(tmp_filename,"_results.txt");
  if(!(results_file = fopen(tmp_filename,"w"))){
  	fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
  	exit(1);
  }
  // Create a file to report results
  strcpy(tmp_filename,basefilename); strcat(tmp_filename,"_tfp.txt");
  if(!(tfp_file = fopen(tmp_filename,"w"))){
  	fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
  	exit(1);
  }
  // Create a file to report minimum pvalues
  strcpy(tmp_filename,basefilename); strcat(tmp_filename,"_minpvals.txt");
  if(!(minpvals_file = fopen(tmp_filename,"w"))){
	fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
	exit(1);
  }

  // Free filename holder
  free(tmp_filename);

  // Additional arguments
  //n_perm = atoi(argv[2]);
  //printf("n_perm %d\n",n_perm);
  //sig_th = atof(argv[3]);
  //printf("sig_th %f\n",sig_th);
  //class_labels_file = argv[4];
  //trans_filename = argv[4];
  //seed_idx = atoi(argv[5]);
  //printf("seed_idx %d\n",seed_idx);
  //int k_ = -1;
  //if(argc > 6)
    //max_eps = atof(argv[6]);
  //printf(stderr, "max_eps %f\n",max_eps);
  //k_ = -1;

  start_measuring_time();

  // Initialise the support of LCM to 1
  LCM_th = 1;
  /* END OF MODIFICATIONS */
  tic = measureTime();
  LCM_problem = LCM_CLOSED;
  LCM_init(input_file);
  toc = measureTime();
  time_LCM_init = toc-tic;

  /* MODIFICATIONS FOR FAST WY ALGORITHIM */
  //Initialise random permutations
  tic = measureTime();
  permutation_init(J,input_file,random_seed);
  toc = measureTime();
  time_permutations = toc-tic;
  if(theta > 0.0){
    LCM_th = (int)(theta * (double)N);
    printf("LCM_th set to %d\n",LCM_th);
  }
  // Initialize MCRapper permutation code
  tic = measureTime();
  mcrapper_init(delta , K, max_eps);
  toc = measureTime();
  time_initialisation_wy = toc-tic;
  /* END OF MODIFICATIONS */
  tic = measureTime();
  LCMclosed();
  toc = measureTime();
  time_threshold_correction = toc-tic;

  // Main part of the code
  LCM_output();
  LCM_end();
  ARY_end(&LCM_Trsact);

  /* MODIFICATION TO KEEP TRACK OF TRANSACTIONS */
  tic = measureTime();
  transaction_keeping_end();
  /* END OF MODIFICATIONS */
  /* MODIFICATIONS FOR FAST WY ALGORITHIM */
  mcrapper_end();
  permutation_end();
  toc = measureTime();
  time_termination = toc-tic;
  // Final time
  t_end = measureTime();
  /* END OF MODIFICATIONS */

  /* MODIFICATIONS FOR CODE PROFILING */
  profileCode();
  /* END OF MODIFICATIONS */

  exit(0);
}


#endif
