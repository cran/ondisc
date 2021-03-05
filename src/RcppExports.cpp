// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// inc_mean_count
void inc_mean_count(NumericVector acc_vect, IntegerVector idxs, IntegerVector umi_counts, double n);
RcppExport SEXP _ondisc_inc_mean_count(SEXP acc_vectSEXP, SEXP idxsSEXP, SEXP umi_countsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type acc_vect(acc_vectSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type umi_counts(umi_countsSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    inc_mean_count(acc_vect, idxs, umi_counts, n);
    return R_NilValue;
END_RCPP
}
// inc_n_entries
void inc_n_entries(IntegerVector acc_vect, IntegerVector idxs);
RcppExport SEXP _ondisc_inc_n_entries(SEXP acc_vectSEXP, SEXP idxsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type acc_vect(acc_vectSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxs(idxsSEXP);
    inc_n_entries(acc_vect, idxs);
    return R_NilValue;
END_RCPP
}
// inc_mean_sq_count
void inc_mean_sq_count(NumericVector acc_vect, IntegerVector idxs, IntegerVector umi_counts, double n);
RcppExport SEXP _ondisc_inc_mean_sq_count(SEXP acc_vectSEXP, SEXP idxsSEXP, SEXP umi_countsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type acc_vect(acc_vectSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type umi_counts(umi_countsSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    inc_mean_sq_count(acc_vect, idxs, umi_counts, n);
    return R_NilValue;
END_RCPP
}
// inc_count
void inc_count(IntegerVector acc_vect, IntegerVector idxs, IntegerVector umi_counts);
RcppExport SEXP _ondisc_inc_count(SEXP acc_vectSEXP, SEXP idxsSEXP, SEXP umi_countsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type acc_vect(acc_vectSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type umi_counts(umi_countsSEXP);
    inc_count(acc_vect, idxs, umi_counts);
    return R_NilValue;
END_RCPP
}
// inc_cell_count_if_feature_condition
void inc_cell_count_if_feature_condition(IntegerVector acc_vect, IntegerVector feature_idxs, IntegerVector cell_idxs, IntegerVector umi_counts, LogicalVector bool_vect);
RcppExport SEXP _ondisc_inc_cell_count_if_feature_condition(SEXP acc_vectSEXP, SEXP feature_idxsSEXP, SEXP cell_idxsSEXP, SEXP umi_countsSEXP, SEXP bool_vectSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type acc_vect(acc_vectSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type feature_idxs(feature_idxsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cell_idxs(cell_idxsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type umi_counts(umi_countsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type bool_vect(bool_vectSEXP);
    inc_cell_count_if_feature_condition(acc_vect, feature_idxs, cell_idxs, umi_counts, bool_vect);
    return R_NilValue;
END_RCPP
}
// index_h5_file
List index_h5_file(const std::string& file_name_in, const std::string& p_name_in, const std::string& idx_name_in, const std::string& umi_counts_name_in, IntegerVector subset_vector, bool logical_mat);
RcppExport SEXP _ondisc_index_h5_file(SEXP file_name_inSEXP, SEXP p_name_inSEXP, SEXP idx_name_inSEXP, SEXP umi_counts_name_inSEXP, SEXP subset_vectorSEXP, SEXP logical_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type file_name_in(file_name_inSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type p_name_in(p_name_inSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type idx_name_in(idx_name_inSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type umi_counts_name_in(umi_counts_name_inSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type subset_vector(subset_vectorSEXP);
    Rcpp::traits::input_parameter< bool >::type logical_mat(logical_matSEXP);
    rcpp_result_gen = Rcpp::wrap(index_h5_file(file_name_in, p_name_in, idx_name_in, umi_counts_name_in, subset_vector, logical_mat));
    return rcpp_result_gen;
END_RCPP
}
// write_data_h5
void write_data_h5(const std::string& file_name_in, const std::string& dataset_name_in, IntegerVector buffer, int start_pos);
RcppExport SEXP _ondisc_write_data_h5(SEXP file_name_inSEXP, SEXP dataset_name_inSEXP, SEXP bufferSEXP, SEXP start_posSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type file_name_in(file_name_inSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type dataset_name_in(dataset_name_inSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type buffer(bufferSEXP);
    Rcpp::traits::input_parameter< int >::type start_pos(start_posSEXP);
    write_data_h5(file_name_in, dataset_name_in, buffer, start_pos);
    return R_NilValue;
END_RCPP
}
// map_memory_to_disk
void map_memory_to_disk(const std::string& file_name_in, IntegerVector m_cell_idxs, const std::string& cell_idxs_name, IntegerVector m_umi_counts, const std::string& umi_counts_name, int n_features, IntegerVector m_row_ptr, IntegerVector f_row_ptr);
RcppExport SEXP _ondisc_map_memory_to_disk(SEXP file_name_inSEXP, SEXP m_cell_idxsSEXP, SEXP cell_idxs_nameSEXP, SEXP m_umi_countsSEXP, SEXP umi_counts_nameSEXP, SEXP n_featuresSEXP, SEXP m_row_ptrSEXP, SEXP f_row_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type file_name_in(file_name_inSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type m_cell_idxs(m_cell_idxsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type cell_idxs_name(cell_idxs_nameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type m_umi_counts(m_umi_countsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type umi_counts_name(umi_counts_nameSEXP);
    Rcpp::traits::input_parameter< int >::type n_features(n_featuresSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type m_row_ptr(m_row_ptrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f_row_ptr(f_row_ptrSEXP);
    map_memory_to_disk(file_name_in, m_cell_idxs, cell_idxs_name, m_umi_counts, umi_counts_name, n_features, m_row_ptr, f_row_ptr);
    return R_NilValue;
END_RCPP
}
// map_memory_to_disk_logical_matrix
void map_memory_to_disk_logical_matrix(const std::string& file_name_in, IntegerVector m_cell_idxs, const std::string& cell_idxs_name, int n_features, IntegerVector m_row_ptr, IntegerVector f_row_ptr);
RcppExport SEXP _ondisc_map_memory_to_disk_logical_matrix(SEXP file_name_inSEXP, SEXP m_cell_idxsSEXP, SEXP cell_idxs_nameSEXP, SEXP n_featuresSEXP, SEXP m_row_ptrSEXP, SEXP f_row_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type file_name_in(file_name_inSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type m_cell_idxs(m_cell_idxsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type cell_idxs_name(cell_idxs_nameSEXP);
    Rcpp::traits::input_parameter< int >::type n_features(n_featuresSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type m_row_ptr(m_row_ptrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f_row_ptr(f_row_ptrSEXP);
    map_memory_to_disk_logical_matrix(file_name_in, m_cell_idxs, cell_idxs_name, n_features, m_row_ptr, f_row_ptr);
    return R_NilValue;
END_RCPP
}
// decrement_idxs
void decrement_idxs(IntegerVector idxs);
RcppExport SEXP _ondisc_decrement_idxs(SEXP idxsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type idxs(idxsSEXP);
    decrement_idxs(idxs);
    return R_NilValue;
END_RCPP
}
// sum_in_place
void sum_in_place(IntegerVector v1, IntegerVector v2);
RcppExport SEXP _ondisc_sum_in_place(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type v2(v2SEXP);
    sum_in_place(v1, v2);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ondisc_inc_mean_count", (DL_FUNC) &_ondisc_inc_mean_count, 4},
    {"_ondisc_inc_n_entries", (DL_FUNC) &_ondisc_inc_n_entries, 2},
    {"_ondisc_inc_mean_sq_count", (DL_FUNC) &_ondisc_inc_mean_sq_count, 4},
    {"_ondisc_inc_count", (DL_FUNC) &_ondisc_inc_count, 3},
    {"_ondisc_inc_cell_count_if_feature_condition", (DL_FUNC) &_ondisc_inc_cell_count_if_feature_condition, 5},
    {"_ondisc_index_h5_file", (DL_FUNC) &_ondisc_index_h5_file, 6},
    {"_ondisc_write_data_h5", (DL_FUNC) &_ondisc_write_data_h5, 4},
    {"_ondisc_map_memory_to_disk", (DL_FUNC) &_ondisc_map_memory_to_disk, 8},
    {"_ondisc_map_memory_to_disk_logical_matrix", (DL_FUNC) &_ondisc_map_memory_to_disk_logical_matrix, 6},
    {"_ondisc_decrement_idxs", (DL_FUNC) &_ondisc_decrement_idxs, 1},
    {"_ondisc_sum_in_place", (DL_FUNC) &_ondisc_sum_in_place, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ondisc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
