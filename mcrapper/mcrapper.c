#ifndef _mcrapper_c_
#define _mcrapper_c_

#define REPORT_INTERVAL 9999999
//#define DEBUG_PERF 1
#define USE_PERM_FILTERING 1
//#define CENTRALIZATION 1
//#define USE_PERM_EARLY_STOP 1

/* LIBRARY INCLUDES */
#include<math.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"
#include"permutation.c"
#include"transaction_keeping.c"
#include"lcm_var.c"
#include"priority_queue.c"
/* CONSTANT DEFINES */

//#define min(a,b) (a < b) ? a : b
//#define max(a,b) (a > b) ? a : b

int min_integer(int a , int b){
  return (a < b) ? a : b;
}

int max_integer(int a , int b){
  return (a > b) ? a : b;
}

double min_double(double a , double b){
  return (a < b) ? a : b;
}

double max_double(double a , double b){
  return (a > b) ? a : b;
}


/* GLOBAL VARIABLES */
FILE* results_file, *minpvals_file, *reports_file, *minsupp_file, *tfp_file;
// Number of samples, N
int N, N_over_2; // consider remove N_over_2
// Number of positive rademacher r.v.
int n; // consider remove
double delta;
double max_pattern_freq;
double max_variance;
double dev_lower_bound;
double approximation_ratio;
double maximum_epsilon;
// maximum deviations for each permutation
double *maxdev;
int *permutations_to_process;
int minimum_deviation_index;
int minimum_deviation;
double avg_rad;
double avg_rad_hyb;
double epsilon_lb;

// store statistics
unsigned long explored_itemsets;
unsigned long patterns_above_topk;
int last_support;
int *supports;
double start_instant;
int last_explored_report;
int tested;

double union_bound_term;
double theta;

// top-k strategy
int K;
heap_t *topk_queue;

// bounds using last tested pattern
int last_precomputation;
int last_processed_support;
int last_processed_as;
int last_processed_max_dev;
int last_processed_min_dev;
int *last_processed_transactions;
int *temp_transaction_list;
double WY_time;
int update_minimum;
// Cell-count counter
int *a_cnt;
int below_delta;
long bound1;
long bound2;
long improved_estimate;
int improved;




void compute_rade_estimate(int k_){
  int j = 0;
  avg_rad = 0.0;
  double min_dev_value = 0.0;
  if(dev_lower_bound > 0.0){
    min_dev_value = (double)LCM_th-1;
  }
  int min_J_k = min_integer((int)J , k_);
  for(j=0; j < min_J_k; j++){
    avg_rad += max_double(maxdev[j] , min_dev_value);
  }
  avg_rad = avg_rad / (double)(min_J_k);
  avg_rad = avg_rad / (double)N;

  // performs centralization
  double dev_labels = 0;
  for(j=0; j < min_J_k; j++){
    dev_labels += 2*n1s[j] - N;
  }
  dev_labels = dev_labels / (double)(min_J_k);
  dev_labels = dev_labels / (double)N;

  avg_rad = avg_rad + dev_labels;
  avg_rad = max_double(avg_rad , 0.0);
}

// Theorem 4.5
void compute_rade_hybrid_estimate(int k_){

  double sample_size = (double)N;
  double perms_ = (double)J;
  double low_freq_bound = sqrt((2.*theta*(log(perms_)+union_bound_term))/sample_size);
  low_freq_bound = low_freq_bound*(double)N;

  int j = 0;
  avg_rad_hyb = 0.0;
  double min_dev_value = low_freq_bound;
  int min_J_k = min_integer((int)J , k_);
  for(j=0; j < min_J_k; j++){
    avg_rad_hyb += max_double(maxdev[j] , min_dev_value);
  }
  avg_rad_hyb = avg_rad_hyb / (double)(min_J_k);
  avg_rad_hyb = avg_rad_hyb / (double)N;

  // performs centralization
  double dev_labels = 0;
  for(j=0; j < min_J_k; j++){
    dev_labels += 2*n1s[j] - N;
  }
  dev_labels = dev_labels / (double)(min_J_k);
  dev_labels = dev_labels / (double)N;

  avg_rad_hyb = avg_rad_hyb + dev_labels;
  avg_rad_hyb = max_double(avg_rad_hyb , 0.0);

}

// Theorem 4.6
double tail_bound_1(double estimated_rademacher){
  double sample_size = (double)N;
  double epsilon_ = 2.*estimated_rademacher + 3.0*sqrt(log(2.0/delta)/(2.0*sample_size));
  return epsilon_;
}

double compute_epsilon_1(){
  compute_rade_estimate(1);
  double epsilon_ = tail_bound_1(avg_rad);
  return epsilon_;
}

// Theorem 3.1
double sbf_tail_bound(double estimated_rademacher){
  double sample_size = (double)N;
  double perms_ = (double)J;
  double z = 0.5;
  double empirical_rademacher = estimated_rademacher + 2.0*z*sqrt(log(4.0/delta)/(2.0*sample_size*perms_));
  double tail_term = sqrt(log(4.0/delta) * (4.0*sample_size*empirical_rademacher + log(4.0/delta)) )/sample_size;
  tail_term += log(4.0/delta)/sample_size;
  tail_term += sqrt(log(4.0/delta)/(2.*sample_size));
  double epsilon_ = 2.*empirical_rademacher + tail_term;
  return epsilon_;
}

double compute_epsilon_sbf(){
  compute_rade_estimate(J);
  double epsilon_ = sbf_tail_bound(avg_rad);
  return epsilon_;
}

double compute_epsilon_sbf_hyb(){
  compute_rade_hybrid_estimate(J);
  double epsilon_ = sbf_tail_bound(avg_rad_hyb);
  return epsilon_;
}

double compute_epsilon_sbf_hyb_1(){
  compute_rade_hybrid_estimate(1);
  double epsilon_ = tail_bound_1(avg_rad_hyb);
  return epsilon_;
}

// Theorem 3.2
double compute_epsilon_sbf_var(){
  compute_rade_estimate(J);
  double sample_size = (double)N;
  double perms_ = (double)J;
  double z = 0.5;
  double log_d = log(4.0/delta);
  double empirical_rademacher = avg_rad + 2.0*z*sqrt(log_d/(2.0*sample_size*perms_));
  empirical_rademacher += (sqrt((4.*sample_size*empirical_rademacher+log_d)*log_d)+log_d)/(2.0*sample_size);
  double epsilon_ = 2.0*empirical_rademacher;
  epsilon_ += sqrt(2.0*log_d * (8.0*empirical_rademacher + max_variance)/sample_size);
  epsilon_ += 2.0*log_d/(3.0*sample_size);
  return epsilon_;
}

double era_tail(double delta , double sample_size){
  double z = 0.5;
  double perms_ = (double)J;
  double log_d = log(4.0/delta);
  return 2.0*z*sqrt(log_d/(2.0*sample_size*perms_));
}

struct bounds_info{
  int frequency;
  //int as;
  int max_as;
  /*int max_as_j;
  int min_as;
  int max_dev;*/
};

void init_bounds(struct bounds_info *current_){
  //printf("Call init_bounds \n");
  current_->frequency = 0;
  //current_->as = 0;
  current_->max_as = 0;
  /*current_->max_as_j = 0;
  current_->min_as = 0;
  current_->max_dev = 0;*/
}

void reset_bounds(struct bounds_info *current_){
  //printf("Call reset_bounds %d , %d\n",current_->frequency,n);
  current_->max_as = current_->frequency;
  /*current_->max_as_j = 0;
  current_->min_as = max(n - (N - current_->frequency) , 0);
  current_->max_dev = current_->frequency;*/
}


/* FUNCTION DECLARATIONS */
int doublecomp(const void*,const void*);

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

void start_measuring_time(){
  start_instant = get_cpu_time();
}

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise MCRapper
 * Input arguments are self-explanatory
 * */
void mcrapper_init(double target_fwer , int k_ , double max_eps_){
  int j; //Loop variable

  delta = target_fwer;
  avg_rad = 0.0;
  avg_rad_hyb = 0.0;

  // Allocate memory for maximum deviations, raising error if it fails
  maxdev = (double *)malloc(J*sizeof(double));
  if(!maxdev){
    fprintf(stderr,"Error in function mcrapper_init: couldn't allocate memory for array maxdev\n");
    exit(1);
  }
  // Array to store the indexes of the permutations to process
  permutations_to_process = (int*)malloc(J*sizeof(int));
  if(!permutations_to_process){
    fprintf(stderr,"Error in function mcrapper_init: couldn't allocate memory for array permutations_to_process\n");
    exit(1);
  }
  // Allocate memory for cell counts, raising an error if it fails
  a_cnt = (int *)malloc(J*sizeof(int));
  if(!a_cnt){
    fprintf(stderr,"Error in function mcrapper_init: couldn't allocate memory for array a_cnt\n");
    exit(1);
  }
  supports = (int *)malloc((N+1) * sizeof(int));
  if(!supports){
    fprintf(stderr,"Error in function mcrapper_init: couldn't allocate memory for array supports\n");
    exit(1);
  }

  // Initialize counts
  for(j=0; j<J; j++){
    a_cnt[j] = 0;
  }
  // Initialise all deviations
  double dev_lower_bound_supp = dev_lower_bound * N;
  for(j=0; j<J; j++){
    if(max_pattern_freq == 1.0){
      maxdev[j] = max_integer(0 , 2*n1s[j] - N );
    }
    else{
      maxdev[j] = 0;
    }
    maxdev[j] = max_integer(maxdev[j] , dev_lower_bound_supp);
  }
  // initialize support counts
  for(j=0; j<=N; j++){
    supports[j] = 0;
  }

  minimum_deviation_index = 0;
  minimum_deviation = 0;
  below_delta = 0;
  explored_itemsets = 0;
  patterns_above_topk = 0;
  last_support = 0;
  last_explored_report = 0;

  K = k_;
  maximum_epsilon = max_eps_;
  epsilon_lb = compute_epsilon_sbf();
  if(maximum_epsilon < 1.0){
    printf("maximum_epsilon %f\n",maximum_epsilon);
    double maximum_rademacher = (maximum_epsilon - epsilon_lb)/2.0;
    printf("maximum value of rademacher %f\n",maximum_rademacher);
    int dev_upper_bound = (int)(maximum_rademacher * (double)N) - 1;
    if(dev_upper_bound > 0 && dev_upper_bound > LCM_th){
      LCM_th = dev_upper_bound;
    }
  }

  reports_file = fopen("reports.txt","w");
  minsupp_file = fopen("minsupp.txt","w");
  last_processed_transactions = (int *)malloc((N+1) * sizeof(int));
  temp_transaction_list =  (int *)malloc((N+1) * sizeof(int));
  topk_queue = (heap_t *)calloc(1, sizeof (heap_t));
  //printf("size of queue: %d\n", topk_queue->len);
  //push(topk_queue, 1.0, 0);
  //pop(topk_queue);
  //printf("size of queue: %d\n", topk_queue->len);
  last_precomputation = -1; last_processed_support = 0;
  last_processed_as = 0; last_processed_max_dev = 0; last_processed_min_dev = 0;
  tested = 0; WY_time = 0.0; bound1 = 0; bound2 = 0; improved_estimate = 0;

  update_minsupp();
  print_debug_report();

}


/* Free all allocated memory and give some output for debugging purposes */
void mcrapper_end(){
  int j, idx_max;

    int min_dev = N;
    int max_dev = 0;
    for (j = 0; j < J; j++) {
        min_dev = min_double(min_dev , maxdev[j]);
        max_dev = max_double(max_dev , maxdev[j]);
    }

    compute_rade_estimate(1);

  // Print results
  fprintf(results_file,"RESULTS\n");
  fprintf(results_file,"\t Estimated Rademacher: %f\n",avg_rad);
  fprintf(results_file,"\t Epsilon lower bound (using standard bounds) %f\n",epsilon_lb);
  fprintf(results_file,"\t Epsilon (n=1) %f\n",compute_epsilon_1());
  fprintf(results_file,"\t Epsilon (no tail bound) %f\n",2.*avg_rad);
  fprintf(results_file,"\t Epsilon (SBF bound) %f\n",compute_epsilon_sbf());
  fprintf(results_file,"\t Epsilon (Var bound) %f\n",compute_epsilon_sbf_var());
  fprintf(results_file,"\t Epsilon (Hyb bound) %f\n",compute_epsilon_sbf_hyb());
  fprintf(results_file,"\t Final LCM support: %d\n",LCM_th);
  fprintf(results_file,"\t Number of explored itemsets: %lu\n",explored_itemsets);
  fprintf(results_file,"\t Number of tested itemsets: %d\n",tested);

  fprintf(results_file,"\t Supports: \n");

  int aj=0;
  for(aj=0; aj < N; aj++){
    fprintf(results_file,"%d : ",aj);
    fprintf(results_file,"%d\n",supports[aj]);
  }

  printf("\n\n");

  double tail_term = era_tail(delta , (double)N);

  printf("Estimated Rademacher: %f\n",avg_rad);
  printf("Estimated hyb Rademacher: %f\n",avg_rad_hyb);
  printf("tail term %f , era %f\n",tail_term,tail_term+avg_rad);
  printf("Min dev: %d , %f\n",min_dev,(double)min_dev / (2. * (double)N));
  printf("Max dev: %d , %f\n",max_dev,(double)max_dev / (2. * (double)N));
  printf("Eps lower bound (using standard bounds) %f\n",epsilon_lb);
  printf("(2*rade estimate with no tail bound) %f\n",2.*avg_rad);
  printf("Epsilon (n=1) %f\n",compute_epsilon_1());
  printf("Epsilon (SBF bound) %f\n",compute_epsilon_sbf());
  printf("Epsilon (Var bound) %f\n",compute_epsilon_sbf_var());
  printf("Epsilon (Hyb bound) %f\n",compute_epsilon_sbf_hyb());
  printf("Epsilon (Hyb bound n=1) %f\n",compute_epsilon_sbf_hyb_1());
  printf("Final LCM support: %d\n",LCM_th);
  printf("Number of explored itemsets: %lu\n",explored_itemsets);
  printf("Number of tested itemsets: %d\n",tested);
  printf("improved_estimate: %ld\n",improved_estimate);
  printf("improved ratio: %f\n",((double)improved_estimate/(double)explored_itemsets));
  printf("bound1: %ld\n",bound1);
  printf("bound2: %ld\n",bound2);
  printf("Runtime for correction (s): %f\n",(get_cpu_time() - start_instant));
  printf("Runtime for WY computations (s): %d\n",(int)WY_time);
  printf("Peak Memory (MB): %d\n",(int)measurePeakMemory());

  // Report the seed for reproducibility
  fprintf(minpvals_file,"\nRANDOM SEED USED FOR KEY GENERATION\n");
  fprintf(minpvals_file,"\t Seed = %u\n",(unsigned)seed);

  // Free allocated memory
  free(maxdev);
  free(permutations_to_process);
  free(a_cnt);
  free(supports);

  free(topk_queue);
  free(last_processed_transactions);
  free(temp_transaction_list);

  // Close files
  fclose(results_file); fclose(minpvals_file);
}

// debug functions printing the deviations
void print_deviations(){
  int j = 0;
  for(j=0; j<J; j++){
    printf("%d,",(int)maxdev[j]);
  }
  printf("\n");
}

// update the frequency threshold for the search
void update_minsupp(){

  minimum_deviation = N;
  int j = 0;
  double new_thr = N;
  double thr_temp = 0.0;
  for(j=0; j<J; j++){
    minimum_deviation = min_double(minimum_deviation , maxdev[j]);
  }
  new_thr = minimum_deviation + 1;

  if(approximation_ratio > 1.){
    new_thr = new_thr * approximation_ratio;
  }

  double epsilon = 0.0;
  if(LCM_th < new_thr){
    LCM_th = new_thr;
    epsilon = compute_epsilon_sbf();
    double tail_term = era_tail(delta , (double)N);
    printf("*  sigma %d, min_dev %d, n1min = %d, est_rade %f, tail_era %f, epsilon %f\n",LCM_th,minimum_deviation,n1s[minimum_deviation_index],avg_rad,tail_term,epsilon);
    printf("   epsilon_lb %f, LCM_th/N %f, est_rade*N %d, eps*N %d\n",epsilon_lb,(double)LCM_th/(double)N,(int)(avg_rad*(double)N),(int)(epsilon*(double)N));
    printf("   explored patterns: %lu, ",explored_itemsets);
    printf("current tested patterns: %d, ",tested);
    printf("tested ratio: %f\n",((double)tested/(double)explored_itemsets));
  }
  if(epsilon > 0.0 && epsilon > maximum_epsilon){
    LCM_th = N+1;
  }

}



int count_common_transactions(int* current_transaction_list , int current_support){
  int common_transactions = 0;
  int i=0; int j=0;
  int cond = 0;
  // count how many transactions are in common between last and current
  while(i<current_support && j<last_processed_support){
    cond = current_transaction_list[i] == last_processed_transactions[j];
    i = (cond) ? i+1 : i;
    j = (cond) ? j+1 : j;
    common_transactions = (cond) ? common_transactions+1 : common_transactions;
    cond = current_transaction_list[i] > last_processed_transactions[j];
    j = (cond) ? j+1 : j;
    i = (cond) ? i : i+1;
  }
  return common_transactions;
}


void print_debug_report(){
  printf("NEW REPORT:\n");
  printf("  sigma %d, est_rade %f, epsilon_lb %f\n",LCM_th,avg_rad,epsilon_lb);
  double epsilon = compute_epsilon_sbf();
  printf("  epsilon %f, LCM_th/N %f, est_rade*N %d, eps*N %d\n",epsilon,(double)LCM_th/(double)N,(int)(avg_rad*(double)N),(int)(epsilon*(double)N));
  printf("  explored patterns: %lu\n",explored_itemsets);
  printf("  current tested patterns: %d\n",tested);
  printf("  tested ratio: %f\n",((double)tested/(double)explored_itemsets));
  printf("  improved_estimate: %ld\n",improved_estimate);
  printf("  improved ratio: %f\n",((double)improved_estimate/(double)explored_itemsets));
  //printf("  last_support:%d\n",last_support);
  printf("  elapsed:%f\n",(get_cpu_time() - start_instant));
  last_explored_report = tested;
  //print_deviations();
}


void process_permutations(struct bounds_info *current_info , int* current_transaction_list , struct bounds_info *parent_info){
  // process the rademacher variables for the current pattern
  //printf("Processing permutations of pattern with support %d\n",current_info->frequency);

  int x = current_info->frequency;
  int i,j;
  char *labels_perm_aux; //Auxiliary pointer
  int curr_max_as = (current_info->max_as < x) ? current_info->max_as : x;//min(current_info->max_as , x);


  #ifdef USE_PERM_FILTERING
  i = 0;
  for(j=0; j<J; j++){
    int current_bound = (curr_max_as < n1s[j]) ? curr_max_as : n1s[j];//min(curr_max_as , n1s[j]);
    if(current_bound > maxdev[j]){
      permutations_to_process[i] = j;
      i++;
    }
  }
  int permutations_to_process_length = i;

  if(permutations_to_process_length > 0){
    if((double)permutations_to_process_length / (double)J < 0.7){
      int current_bound;
      int perm_idx;
      #ifdef DEBUG_PERF
      printf("processing labels\n");
      #endif
      for(i=0; i<x; i++){
        #ifdef DEBUG_PERF
        printf("    current_transaction_list[i] %d",current_transaction_list[i]);
        #endif
        labels_perm_aux = labels_perm[current_transaction_list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
        for(j=0; j<permutations_to_process_length; j++){
          perm_idx = permutations_to_process[j];
          a_cnt[perm_idx] += labels_perm_aux[perm_idx];

          #ifdef USE_PERM_EARLY_STOP
          if((i+1)%10==0){ // don't check too often
            current_bound = 2*(a_cnt[permutations_to_process[j]] + min_integer(curr_max_as,x - (i+1)) ) - x;
            if(current_bound <= maxdev[j]){
              // remove current permutation from further consideration
              permutations_to_process[j] = permutations_to_process[permutations_to_process_length-1];
              permutations_to_process_length--;
            }
          }
          #endif

        }
      }
    }
    else{
      for(i=0; i<x; i++){
        labels_perm_aux = labels_perm[current_transaction_list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
        for(j=0; j<J; j++){
          a_cnt[j] += labels_perm_aux[j];
        }
      }
    }
  }
  #endif
  #ifndef USE_PERM_FILTERING
  // standard loop to update values of a_cnt
  for(i=0; i<x; i++){
    labels_perm_aux = labels_perm[current_transaction_list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
    for(j=0; j<J; j++){
      a_cnt[j] += labels_perm_aux[j];
    }
  }
  #endif

  // the current will be the new last processed pattern
  #ifdef ENABLE_LAST_TESTED_BOUNDS
  for(i=0; i<x; i++) {last_processed_transactions[i] = current_transaction_list[i]; }
  #endif
  last_processed_support = x;
  last_processed_min_dev = x;
  last_processed_max_dev = 0;

  update_minimum = 0;
  improved = 0;
  int any_above = 0;
  int max_dev_above = 0;
  int current_deviation = 0;
  int current_dev_compl = 0;
  int max_j_index = 0;
  int max_dev = N;
  int current_bound;
  curr_max_as = 0;

  for(j=0; j<J; j++){
    current_deviation = 2 * a_cnt[j] - x;
    curr_max_as = (curr_max_as < a_cnt[j]) ? a_cnt[j] : curr_max_as;
    improved += (current_deviation > maxdev[j]);
    maxdev[j] = (current_deviation > maxdev[j]) ? current_deviation : maxdev[j];
    a_cnt[j] = 0;
  }
  update_minimum = (improved > 0);
  improved_estimate += (improved > 0);
  current_info->max_as = curr_max_as;

}

#define ENABLE_BOUNDS 1
//#define ENABLE_LAST_TESTED_BOUNDS 1

int current_pattern_is_testworthy(struct bounds_info* parent_info , struct bounds_info* current_info ,  int* current_transaction_list){

  // don't use bounds, always process the pattern
  #ifndef ENABLE_BOUNDS
  return 0;
  #endif

  int to_ret = 0;
  int j;
  int any_above = 0;
  int x = current_info->frequency;
  int current_max_as = (current_info->max_as < x) ? current_info->max_as : x;
  int max_as_j_;
  int upp_bound_dev_j = 0;
  int max_bound_dev_j = 0;
  #ifdef DEBUG_PERF
  fflush(stdout);
  printf("untest check: x %d , max_as %d\n",x,current_max_as);
  #endif
  for(j=0; j < J; j++){
    max_as_j_ = (current_max_as < n1s[j]) ? current_max_as : n1s[j];
    upp_bound_dev_j = 2*max_as_j_ - x;
    any_above += ( upp_bound_dev_j > maxdev[j] );
    #ifdef DEBUG_PERF
    max_bound_dev_j = max_integer(max_bound_dev_j,upp_bound_dev_j);
    printf("   j %d , maxdev[j] %f , n1[j] %d , max_as_j %d , upp_bound_dev_j %d \n",j,maxdev[j],n1s[j],max_as_j_,upp_bound_dev_j);
    #endif
  }
  if(any_above == 0){
    to_ret = 1;
    bound2++;
  }
  #ifdef DEBUG_PERF
  printf("   any_above %d , max_bound_dev_j %d , toret %d\n",any_above,max_bound_dev_j,to_ret);
  fflush(stdout);
  #endif

  // compute a bound on the current pattern using the last processed pattern
  // possibly broken , to be fixed
  #ifdef ENABLE_LAST_TESTED_BOUNDS
  //if (current_transaction_list && ( ( last_processed_max_dev + (current_info->frequency - count_common_transactions(current_transaction_list , current_info->frequency))) < LCM_th ) )
  //  return 1;
  #endif

  return to_ret;

}

/* Process a solution involving the bitmap represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
void bm_process_solution(int x, int item, int *mask , struct bounds_info* parent_info , struct bounds_info* current_info){
  //printf("bm_process_solution: supp %d\n",x);
  //printf("  current tested patterns:%d\n",tested);
  int i,j;//Loop iterators

  if(max_pattern_freq == 0.0){
    fprintf(tfp_file,"%d - %d\n",item,x);
  }

  //double temp_time_ = get_cpu_time();

  if(parent_info->frequency <= 0){
    parent_info->frequency = N;
    parent_info->max_as = N;
  }
  current_info->max_as = (x <= parent_info->max_as) ? x : parent_info->max_as;//min(x , parent_info->max_as);
  current_info->frequency = x;

  if(x != current_info->frequency){
    printf("problem in bm_process_solution: supp %d, current supp %d\n",x,current_info->frequency);
  }
  if(x < N && x >= parent_info->frequency){
    printf("problem in bm_process_solution: supp %d, parent supp %d\n",x,parent_info->frequency);
  }

  ++explored_itemsets;
  ++supports[x];

  last_support = x;

  // Sanity-check
  if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);

  if( ((double)x)/((double)N) > max_pattern_freq) return;

  if( current_pattern_is_testworthy(parent_info , current_info , current_trans.list) ) return;

  #ifdef PROFILE_MINING_OLD
  ticp = measureTime();
  #endif


  ++tested;

  process_permutations(current_info , current_trans.list , parent_info);

  /* Finally, check if the minimum constraint is still satisfied, if not compute it again */
  if(update_minimum == 1) {
    update_minsupp();
    // Correct possible corruption of LCM data structures due to unexpected change in minimum support
    for(i=0; i<item; i++){
      //printf("Item %d, Frq %d, Current_th %d\n",i,LCM_Ofrq[i],LCM_th);
      if(LCM_Ofrq[i]==(LCM_th-1)){
        //printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",i,item,LCM_th);
        LCM_BM_occurrence_delete(i);
        *mask &= ~BITMASK_1[i];
        //printf("Problem fixed!\n");
      }
    }
  }
  #ifdef PROFILE_MINING_OLD
  tocp = measureTime(); time_minpval += tocp-ticp;
  #endif

  if(tested - last_explored_report > REPORT_INTERVAL){
    print_debug_report();
  }


}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x , struct bounds_info* parent_info , struct bounds_info* current_info){
  //printf("process_solution0: supp %d\n",x);
  //printf("  current tested patterns:%d\n",tested);
  int i,j;//Loop iterators

  if(max_pattern_freq == 0.0){
    fprintf(tfp_file,"%d - %d\n",0,x);
  }

  if(parent_info->frequency <= 0){
    parent_info->frequency = N;
    parent_info->max_as = N;
  }
  current_info->max_as = (x <= parent_info->max_as) ? x : parent_info->max_as;//min(x , parent_info->max_as);
  current_info->frequency = x;

  if(x != current_info->frequency){
    printf("problem in bm_process_solution: supp %d, current supp %d\n",x,current_info->frequency);
  }
  if(x < N && x >= parent_info->frequency){
    printf("problem in bm_process_solution: supp %d, parent supp %d\n",x,parent_info->frequency);
  }

  ++explored_itemsets;
  ++supports[x];

  last_support = x;

  // Sanity-check
  if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);

  if( ((double)x)/((double)N) > max_pattern_freq) return;

  if( current_pattern_is_testworthy(parent_info , current_info , bm_trans_list[1].list) ) return;

  #ifdef PROFILE_MINING_OLD
  ticp = measureTime();
  #endif

  ++tested;

  process_permutations(current_info , bm_trans_list[1].list , parent_info);
  /* Finally, check if the minimum constraint is still satisfied, if not compute it again */
  if(update_minimum == 1) {
    update_minsupp();
  }

  #ifdef PROFILE_MINING_OLD
  tocp = measureTime(); time_minpval += tocp-ticp;
  #endif

  if(tested - last_explored_report > REPORT_INTERVAL){
    print_debug_report();
  }

}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution(int x, TRANS_LIST *L, int item, int *mask , struct bounds_info* parent_info , struct bounds_info* current_info){
  //printf("ary_process_solution: supp %d\n",x);
  //printf("  current tested patterns:%d\n",tested);
  int i,j;//Loop iterator
  int aux; //Auxiliary counter
  int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list

  //double temp_time_ = get_cpu_time();

  if(max_pattern_freq == 0.0){
    fprintf(tfp_file,"%d - %d\n",item,x);
  }

  if(parent_info->frequency <= 0){
    parent_info->frequency = N;
    parent_info->max_as = N;
  }
  current_info->max_as = (x <= parent_info->max_as) ? x : parent_info->max_as;//min(x , parent_info->max_as);
  current_info->frequency = x;

  if(x != current_info->frequency){
    printf("problem in bm_process_solution: supp %d, current supp %d\n",x,current_info->frequency);
  }
  if(x < N && x >= parent_info->frequency){
    printf("problem in bm_process_solution: supp %d, parent supp %d\n",x,parent_info->frequency);
  }

  ++explored_itemsets;
  ++supports[x];

  last_support = x;


  if( current_pattern_is_testworthy(parent_info , current_info , NULL) ) return;

  #ifdef PROFILE_MINING_OLD
  ticp = measureTime();
  #endif

  // Sanity-check (this one is more complicated due to the way the transactions are stored)
  aux = 0;
  for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
    end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
    for(ptr = L->ptr[*t];ptr < end_ptr;ptr++) aux++;
  }
  if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);


      j=0;
      for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
        end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
        for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
          temp_transaction_list[j] = *ptr;
          ++j;
        }
      }

    if( ((double)x)/((double)N) > max_pattern_freq) return;

    #ifdef ENABLE_LAST_TESTED_BOUNDS
    if( current_pattern_is_testworthy(parent_info , current_info , temp_transaction_list) ) return;
    #endif

    ++tested;

    process_permutations(current_info , temp_transaction_list , parent_info);
    /* Finally, check if the minimum constraint is still satisfied, if not compute it again */
    if(update_minimum == 1) {
      update_minsupp();
      // Correct possible corruption of LCM data structures due to unexpected change in minimum support
      for(j=0; j<LCM_BM_MAXITEM; j++){
        //printf("Item %d, Frq %d, Current_th %d\n",j,LCM_Ofrq[j],LCM_th);
        if(LCM_Ofrq[j]==(LCM_th-1)){
          //printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",j,item,LCM_th);
          LCM_BM_occurrence_delete(j);
          *mask &= ~BITMASK_1[j];
          //printf("Problem fixed!\n");
        }
      }
    }

  #ifdef PROFILE_MINING_OLD
  tocp = measureTime(); time_minpval += tocp-ticp;
  #endif

  if(tested - last_explored_report > REPORT_INTERVAL){
    print_debug_report();
  }


}

/* AUXILIARY FUNCTIONS */
// Comparison function used by quicksort implementation in C
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
}

#endif
