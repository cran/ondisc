########
# Test 1
########
test_that("extract by gene id", {
  for (i in 1:n_datasets) {
    # Load data
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    # Set row and column names for Matrix object
    row.names(Mat) <- get_feature_ids(on_disc_mat)
    colnames(Mat) <- get_cell_barcodes(on_disc_mat)
    for (j in 1:n_reps) {
      # Test arbitrary subsets of rows and columns
      subset_size_col <- sample(1:(ceiling(ncol(Mat)/30)), 1)
      subset_size_row <- sample(1:(ceiling(nrow(Mat)/30)), 1)
      col_names <- paste0("cell_", sample(x = 1:ncol(Mat), size = subset_size_col))
      row_names <- paste0("ENSG000", sample(x = 1:nrow(Mat), size = subset_size_row))
      compare_Mat_on_disc_extract(Mat, on_disc_mat, col_names, row_names)
    }
  }
})

########
# Test 2
########
test_that("extract by logical vector", {
  for (i in 1:n_datasets) {
    # Load data
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    for (j in 1:n_reps) {
      col_subset <- sample(x = c(FALSE, TRUE), size = ncol(Mat), replace = TRUE, prob = c(0.98, 0.02))
      row_subset <- sample(x = c(FALSE, TRUE), size = nrow(Mat), replace = TRUE, prob = c(0.98, 0.02))
      if (any(col_subset) && any(row_subset)) {
        compare_Mat_on_disc_extract(Mat = Mat, on_disc_mat = on_disc_mat, col_idxs = col_subset, row_idxs = row_subset)
      }
    }
  }
})

########
# Test 3
########
test_that("extract by gene id on subset matrix", {
  for (i in 1:n_datasets) {
    # Load data
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    row.names(Mat) <- paste0("ENSG000", 1:nrow(Mat))
    colnames(Mat) <- paste0("cell_", 1:ncol(Mat))
    for (j in 1:n_reps) {
      row_names <- sample(get_feature_ids(on_disc_mat), sample(1:nrow(on_disc_mat), 1))
      col_names <- sample(get_cell_barcodes(on_disc_mat), sample(1:ncol(on_disc_mat), 1))
      # Do the subsets
      Mat_col_sub <- Mat[,col_names,drop=FALSE]
      Mat_row_sub <- Mat[row_names,,drop=FALSE]
      Mat_sub <- Mat[row_names, col_names,drop=FALSE]
      on_disc_col_sub <- on_disc_mat[,col_names]
      on_disc_row_sub <- on_disc_mat[row_names,]
      on_disc_sub <- on_disc_mat[row_names,col_names]
      # Obtain the extraction subsets
      col_names_extract <- sample(col_names, sample(1:ceiling(length(col_names)/20), 1))
      row_names_extract <- sample(row_names, sample(1:ceiling(length(row_names)/20), 1))
      # Test the extracts
      compare_Mat_on_disc_extract(Mat_col_sub, on_disc_col_sub, col_names_extract, row_names)
      compare_Mat_on_disc_extract(Mat_row_sub, on_disc_row_sub, col_names, row_names_extract)
      compare_Mat_on_disc_extract(Mat_sub, on_disc_sub, col_names_extract, row_names_extract)
    }
  }
})

########
# Test 4
########
test_that("Extract by logical vector on subset matrix", {
  for (i in 1:n_datasets) {
    # Load data
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    for (j in 1:n_reps) {
      col_idxs <- sample(1:ncol(on_disc_mat), sample(1:ceiling(ncol(on_disc_mat)/30), 1))
      row_idxs <- sample(1:nrow(on_disc_mat), sample(1:ceiling(nrow(on_disc_mat)/30), 1))
      Mat_sub <- Mat[1:nrow(Mat) %in% row_idxs, 1:ncol(Mat) %in% col_idxs, drop = FALSE]
      on_disc_sub <- on_disc_mat[1:nrow(on_disc_mat) %in% row_idxs, 1:ncol(on_disc_mat) %in% col_idxs]
      col_sub_idx <- sample(x = c(TRUE, FALSE), size = ncol(Mat_sub), replace = TRUE)
      row_sub_idx <- sample(x = c(TRUE, FALSE), size = nrow(Mat_sub), replace = TRUE)
      if (any(row_sub_idx) && any(col_sub_idx)) {
        compare_Mat_on_disc_extract(Mat = Mat_sub, on_disc_mat = on_disc_sub, col_idxs = col_sub_idx, row_idxs = row_sub_idx)
      }
    }
  }
})

########
# Test 5
########
# Test illegal requests
test_that("subset by gene id", {
  for (i in 1:n_datasets) {
    # Load data
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    # Elements not present
    expect_error(on_disc_mat[,c("cell_1", "bogus name")])
    expect_error(on_disc_mat[c("ENSG0001", "bogus name"),])
    # Vector too long
    expect_error(on_disc_mat[, c(rep(TRUE, ncol(on_disc_mat)), TRUE)])
    expect_error(on_disc_mat[c(rep(TRUE, nrow(on_disc_mat)), TRUE),])
    # No cell or genes in on_dist_matrix
    expect_error(on_disc_mat[,FALSE][[1,]])
    expect_error(on_disc_mat[FALSE,][[,1]])
  }
})
