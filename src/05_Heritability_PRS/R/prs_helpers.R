# Checks specific to the original chromosome 22 PRS exercise.
prs_require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(x)
}

prs_check_ids <- function(ids, label) {
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(label, " must contain nonmissing, nonempty, unique IDs.")
  }
  invisible(ids)
}

prs_as_numeric <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(value))) {
    stop(label, " contains missing, nonnumeric, or nonfinite values.")
  }
  value
}

prs_read_fam <- function(path) {
  fam <- data.table::fread(path, header = FALSE, colClasses = "character")
  if (ncol(fam) != 6L) stop("The target .fam file must have six columns.")
  data.table::setnames(fam, c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
  prs_check_ids(fam$IID, "Target .fam IID")
  if (anyNA(fam$FID) || any(!nzchar(fam$FID))) {
    stop("Target .fam FID contains missing values.")
  }
  fam
}

prs_read_outcomes <- function(path, fam) {
  outcome <- data.table::fread(path, colClasses = "character")
  prs_require_columns(outcome, c("ID", "y"), "y_out")
  prs_check_ids(outcome$ID, "y_out ID")
  if (nrow(outcome) != 20000L) {
    stop("The original exercise requires exactly 20,000 rows in y_out; found ",
         nrow(outcome), ".")
  }
  if (!setequal(outcome$ID, fam$IID)) {
    stop("y_out ID must match the target .fam IID set exactly. ",
         "Do not replace this check with a row-order join.")
  }
  outcome[, y := prs_as_numeric(y, "y_out y")]
  outcome[, original_row := .I]
  outcome[, subset := rep(c("tuning", "validation"), each = 10000L)]
  for (group in c("tuning", "validation")) {
    if (stats::var(outcome[subset == group, y]) == 0) {
      stop("The ", group, " outcome is constant; R-squared is undefined.")
    }
  }
  outcome
}

prs_read_bim <- function(path) {
  variants <- data.table::fread(path, header = FALSE, colClasses = "character")
  if (ncol(variants) != 6L) stop("The target .bim file must have six columns.")
  data.table::setnames(variants, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  prs_check_ids(variants$SNP, "Target .bim SNP")
  if (anyNA(variants$A1) || anyNA(variants$A2) ||
      any(!nzchar(variants$A1)) || any(!nzchar(variants$A2)) ||
      any(variants$A1 == variants$A2)) {
    stop("Target .bim alleles must be nonmissing and distinct at each variant.")
  }
  variants
}

prs_align_scores <- function(path, outcome, fam) {
  score <- data.table::fread(path, colClasses = "character")
  # PLINK 2 prefixes the first ID column with '#'.
  data.table::setnames(score, sub("^#", "", names(score)))
  prs_require_columns(score, c("IID", "SCORE1_SUM"), "PLINK 2 score file")
  prs_check_ids(score$IID, "PLINK 2 score IID")
  if (!setequal(score$IID, outcome$ID)) {
    stop("Score IID and y_out ID sets differ in ", path, ".")
  }
  if ("FID" %in% names(score)) {
    expected_fid <- fam$FID[match(score$IID, fam$IID)]
    if (anyNA(score$FID) || any(score$FID != expected_fid)) {
      stop("Score FID/IID combinations do not match the target .fam file.")
    }
  }
  score[, SCORE1_SUM := prs_as_numeric(SCORE1_SUM, "SCORE1_SUM")]
  aligned <- data.table::copy(outcome)
  aligned[, score := score$SCORE1_SUM[match(ID, score$IID)]]
  aligned
}
