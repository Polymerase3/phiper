# Generates delta_regression.rds for test-delta_regression.R.
# Run from the package root: Rscript tests/testthat/fixtures/make_delta_regression.R

set.seed(1)
n_subj <- 300
n_pep <- 50
subjects <- sprintf("s%03d", seq_len(n_subj))
peptides <- sprintf("p%02d", seq_len(n_pep))
group <- rep(c("A", "B"), each = n_subj / 2)

# prevalence 2-40%; every second peptide is shifted up in group B and at T2
base_p <- seq(0.02, 0.40, length.out = n_pep)
shift <- rep(c(0, 0.05), length.out = n_pep)

data <- expand.grid(
  peptide_id = peptides,
  subject_id = subjects,
  timepoint = c("T1", "T2"),
  stringsAsFactors = FALSE
)
pi <- match(data$peptide_id, peptides)
si <- match(data$subject_id, subjects)
p <- base_p[pi] + shift[pi] * (group[si] == "B") + shift[pi] * (data$timepoint == "T2")
data$exist <- rbinom(nrow(data), 1, p)
data$group <- group[si]
data$sample_id <- paste(data$subject_id, data$timepoint, sep = "_")

peptide_library <- data.frame(
  peptide_id = peptides,
  species = rep(sprintf("species_%d", 1:5), each = 10)
)

saveRDS(
  list(data = data, peptide_library = peptide_library),
  "tests/testthat/fixtures/delta_regression.rds"
)
