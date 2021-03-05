#include <Rcpp.h>
using namespace Rcpp;

//' @title increment mean count
//' @param acc_vect a numeric accumulator vector
//' @param idxs 0-based integer vector of indexes (e.g., feature indexes)
//' @param umi_counts integer vector of UMI counts
//' @param n a double indicating the total number of elements to divide by (e.g., number of cells)
//' @noRd
// [[Rcpp::export]]
void inc_mean_count(NumericVector acc_vect, IntegerVector idxs, IntegerVector umi_counts, double n) {
  for (int i = 0; i < idxs.size(); i++) {
    acc_vect[idxs[i]] += umi_counts[i]/n;
  }
}


//' @title increment n entries
//' @param acc_vect an integer accumulator vector
//' @param idxs 0-based integer vector of indexes
//' @noRd
// [[Rcpp::export]]
void inc_n_entries(IntegerVector acc_vect, IntegerVector idxs) {
  for (int i = 0; i < idxs.size(); i ++) {
    acc_vect[idxs[i]] ++;
  }
}


//' @title increment mean squared count
//' @param acc_vect a numeric accumulator vector
//' @param idxs 0-based integer vector of indexes (e.g., feature indexes)
//' @param umi_counts integer vector of UMI counts
//' @param n a double indicating the total number of elements to divide by (e.g., number of cells)
//' @noRd
// [[Rcpp::export]]
void inc_mean_sq_count(NumericVector acc_vect, IntegerVector idxs, IntegerVector umi_counts, double n) {
  int curr_count;
  for (int i = 0; i < idxs.size(); i++) {
    curr_count = umi_counts[i];
    acc_vect[idxs[i]] += (curr_count * curr_count)/n;
  }
}


//' @title increment count
//' @param acc_vect an integer accumulator vector
//' @param idxs 0-based integer vector of indexes
//' @param umi_counts integer vector of UMI counts
//' @noRd
// [[Rcpp::export]]
void inc_count(IntegerVector acc_vect, IntegerVector idxs, IntegerVector umi_counts) {
  for (int i = 0; i < idxs.size(); i++) {
    acc_vect[idxs[i]] += umi_counts[i];
  }
}


//' @title increment count (if condition on feature holds)
//' @param acc_vect an integer accumulator vector
//' @param feature_idxs 0-based integer vector of feature indexes
//' @param cell_idxs 0-based integer vector of cell_idxs
//' @param umi_counts integer vector of UMI counts
//' @param bool_vect a logical vector indicating whether to increment acc_vect; should be same length as n_features
//' @noRd
// [[Rcpp::export]]
void inc_cell_count_if_feature_condition(IntegerVector acc_vect, IntegerVector feature_idxs, IntegerVector cell_idxs, IntegerVector umi_counts, LogicalVector bool_vect) {
  for (int i = 0; i < cell_idxs.size(); i++) {
    if (bool_vect[feature_idxs[i]]) {
      acc_vect[cell_idxs[i]] += umi_counts[i];
    }
  }
}
