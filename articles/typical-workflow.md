# Tutorial: Typical phiper workflow

## General idea

This tutorial is supposed to provide an exemplary workflow for a
cross-sectional analysis of unpaired PhIP-Seq data using the `phiper`
package. It’s supposed to serve as a railing along `phiper`’s individual
modules to help you get started. It’s not meant as a comprehensive
collection of all features.

For instructions to install `phiper`, please consult the
[documentation](https://polymerase3.github.io/phiper/index.html).

`phiper`’s modules are explained in detail on their own documentation
pages, linked here where appropriate. All of `phiper`’s functions also
come with their own help page, so please call `?function_name` for
focused usage instructions.

We provide a dummy dataset for you to play with and we use it here to
showcase `phiper`’s analysis and plotting functions.

## Load packages

You will want to start by loading the necessary packages. `phiperio`
handles input/output of phiper-related data, while `phiper` contains
functions for analysis and plotting. We also use `dplyr` and `ggplot2`
here for quality of life.

``` r

library(phiperio)
library(phiper)
library(dplyr)
library(ggplot2)
```

## Load data

Next, you will need to load your data. There are different ways to do
this. For details, please consult `phiperio`’s [Github
repository](https://github.com/Polymerase3/phiperio). For this tutorial
we assume that you have received a set of `.csv` files (one for each
sample), containing a list of enriched peptides.

Note that peptides which didn’t show enrichment within a sample will not
be listed in that sample’s `.csv` file.

Please save the entire set of `.csv` files that you want to load into
one folder. Then `phiperio`’s
[`convert_standard()`](https://polymerase3.github.io/phiperio/reference/convert_standard.html)
can batch-import them into a `phip-data` object.
[`convert_standard()`](https://polymerase3.github.io/phiperio/reference/convert_standard.html)
can also be used to load a `.parquet` file.

For this tutorial, we have prepared a dummy dataset bundled with the
package. Load it like so:

``` r

vig_dir <- system.file("extdata/typical-workflow", package = "phiper")

pd <- convert_standard(
  data_long_path = file.path(vig_dir, "enrichment_files"),
  peptide_id = "ID",                   # column in the .csv files for peptide IDs
  sample_id_from_filenames = TRUE,     # use file names as sample IDs
  auto_expand = TRUE                   # expand to include all library peptides
)
#> [19:51:22] INFO  Constructing <phip_data> object
#>                  -> create_data()
#> [19:51:22] INFO  Fetching peptide metadata library via get_peptide_library()
#> [19:51:22] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(force_refresh = FALSE)
#> [19:51:22] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - table: peptide_meta
#> [19:51:22] OK    Using cached download (SHA-256 match)
#> [19:51:25] OK    Download complete and loaded into R
#> [19:51:30] INFO  Importing sanitized metadata into DuckDB cache...
#> [19:51:32] OK    peptide_meta table created in DuckDB cache
#> [19:51:32] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 9.359s
#> [19:51:32] OK    Peptide metadata acquired
#> [19:51:32] INFO  Validating <phip_data>
#>                  -> validate_phip_data()
#> [19:51:32] INFO  Checking structural requirements (shape & mandatory columns)
#> [19:51:32] INFO  Checking outcome family availability (exist / fold_change /
#>                  raw_counts)
#> [19:51:32] INFO  Checking collisions with reserved names
#>                    - subject_id, sample_id, timepoint, peptide_id, exist,
#>                      fold_change, counts_input, counts_hit
#> [19:51:32] INFO  Ensuring all columns are atomic (no list-cols)
#> [19:51:32] INFO  Checking key uniqueness
#> [19:51:32] INFO  Validating value ranges & types for outcomes
#> [19:51:32] INFO  Assessing sparsity (NA/zero prevalence vs threshold)
#>                    - warn threshold: 50%
#> [19:51:32] INFO  Checking peptide_id coverage against peptide_library
#> [19:51:32] INFO  Checking full grid completeness (peptide * sample)
#> [19:51:32] INFO  Auto-expanding to full grid via expand_data()
#>                    - add_exist = TRUE
#>                    - exist_col = "exist"
#> [19:51:32] INFO  Expanding <phip_data> to full grid
#>                  -> updating x$data_long
#> [19:51:32] INFO  Expanding to full key * id grid
#>                  -> keys: 'sample_id'; id: 'peptide_id'
#> [19:51:32] INFO  Type probe on lazy table
#>                  -> collect(head 0)
#> [19:51:32] INFO  Building Cartesian product of keys and ids
#> [19:51:32] INFO  Detecting per-key constant (recyclable) columns
#>                    - candidates: fold_change, neglogp, padj, input, count
#> [19:51:32] OK    Column split decided
#>                    - recyclable: <none>
#>                    - non-recyclable: fold_change, neglogp, padj, input, count
#> [19:51:32] INFO  Preparing fill defaults for introduced rows
#>                    - numeric/integer: fold_change, neglogp, padj, input, count
#>                    - logical: <none>
#> [19:51:32] INFO  Applying user-provided fill overrides
#>                    - overrides: exist, fold_change, input_count, hit_count,
#>                      counts_input, counts_hit
#> [19:51:32] INFO  Adding existence flag column
#>                    - column: "exist"
#> [19:51:32] OK    Expanding to full key * id grid - done
#>                  -> elapsed: 0.156s
#> [19:51:33] INFO  Registering expanded table back to DB
#>                    - name: 'data_long'
#>                    - materialise_table: TRUE
#> [19:51:33] INFO  Registering lazy table
#>                  -> name: 'data_long'; as TABLE
#> [19:51:33] INFO  Materialising via dplyr::compute()
#> [19:51:33] OK    Registering lazy table - done
#>                  -> elapsed: 0.251s
#> [19:51:33] OK    Expanding <phip_data> to full grid - done
#>                  -> elapsed: 1.225s
#> [19:51:33] OK    Auto-expansion complete; grid is now full
#> [19:51:33] OK    Validating <phip_data> - done
#>                  -> elapsed: 1.604s
#> [19:51:33] OK    Constructing <phip_data> object - done
#>                  -> elapsed: 10.966s
```

Since we are dealing with extremely large tables, `phip_data` objects
rely on DuckDB, so the data is not stored in memory but on disk instead.
Calls and manipulations of this object will be translated into database
queries. This can be a bit awkward at times, but you can basically do
the same types of operations as on regular R objects (e.g.,
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html),
[`select()`](https://dplyr.tidyverse.org/reference/select.html), `%>%`,
etc.). This is how calling this object looks like:

``` r

pd
#> ── <phip_data> ───────────────────────────────────────────────────────────────── 
#> 
#> counts (first 5 rows): 
#> # A tibble: 5 × 8
#>   sample_id peptide_id     fold_change neglogp  padj input count exist
#>   <chr>     <chr>                <dbl>   <dbl> <dbl> <dbl> <dbl> <int>
#> 1 B_20      agilent_196462     13449.     8.20 0.002   0.1    18     1
#> 2 B_20      agilent_165933        97.7    7.02 0.029  13      17     1
#> 3 B_20      twist_83556           61.7    7.35 0.013  23      19     1
#> 4 B_20      twist_75961         1644.     9.71 0       1      22     1
#> 5 B_20      twist_54336           88.7    7.66 0.007  16      19     1
#> 
#> table size: 89,000 rows x 8 columns
#> 
#> peptide library preview (first 5 rows): 
#> # A tibble: 5 × 8
#>   peptide_id Fullname                    species genus family order class common
#>   <chr>      <chr>                       <chr>   <chr> <chr>  <chr> <chr> <chr> 
#> 1 agilent_1  Chromodomain-helicase-DNA-… Homo s… Homo  Homin… Prim… Mamm… Human 
#> 2 agilent_2  integral membrane protein   Mycoba… Myco… Mycob… Myco… Acti… NA    
#> 3 agilent_3  hypothetical protein (6/16… Mycoba… Myco… Mycob… Myco… Acti… NA    
#> 4 agilent_4  envelope protein (5/8) & a… Orthof… Orth… Flavi… Amar… Flas… JEV   
#> 5 agilent_5  Myosin-7 & beta-myosin hea… Homo s… Homo  Homin… Prim… Mamm… Human 
#> ... plus 36 more columns
#> 
#> library size: 357,190 rows x 44 columns
#> 
#> meta flags: 
#>   con:            <duckdb_connection>
#>   longitudinal:   FALSE
#>   exist:          TRUE
#>   fold_change:    TRUE
#>   raw_counts:     FALSE
#>   extra_cols:     neglogp, padj, input, count
#>   peptide_con:    <duckdb_connection>
#>   materialise_table: TRUE
#>   finalizer_env:  <environment>
#>   full_cross:     TRUE
#>   exist_prop:     0.230280898876404
```

You can extract a regular tibble from a `phip_data` object using
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html), but
larger datasets will have to be filtered down first, as they won’t fit
into memory.

`phip_data` objects also come with annotation information of all
peptides in the library.
[`get_peptide_library()`](https://polymerase3.github.io/phiperio/reference/get_peptide_library.html)
will return a `phip_data` object with these annotations. You may want to
turn that into a tibble for convenient exploration:

``` r

peplib <- get_peptide_library() %>%
  collect()
#> [19:51:33] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(force_refresh = FALSE)
#> [19:51:34] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - table: peptide_meta
#> [19:51:34] OK    Using cached download (SHA-256 match)
#> [19:51:35] OK    Download complete and loaded into R
#> [19:51:40] INFO  Importing sanitized metadata into DuckDB cache...
#> [19:51:42] OK    peptide_meta table created in DuckDB cache
#> [19:51:42] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 8.249s
```

Let’s now load in the metadata for our samples and include them into our
`pd` object. The metadata for our dummy dataset are entirely fictional.
This metadata table contains a column called `treatment` assigning the
samples into either `treat_A` or `treat_B`. We will be using this column
for comparisons throughout this tutorial. Across different functions,
comparisons are realized by specifying the column name in the
`group_cols` argument (i.e., `group_cols = "treatment"`).

``` r

metadata <- read.delim(file.path(vig_dir, "metadata.txt")) %>%
  tibble()

pd_with_metadata <- left_join(pd, metadata, by = "sample_id", copy = TRUE) # copy = TRUE is necessary because of DuckDB
```

`pd_with_metadata` will be the main input for most of `phiper`’s
functions in this tutorial.

## Alpha diversity

Please refer to the vignette [Alpha Diversity Analysis in
PhIP-seq](https://polymerase3.github.io/phiper/articles/alpha-diversity.html)
for detailed information.

As a first overview, let’s plot the distribution of enriched peptide
counts across our samples.

``` r

p_enrichment_counts <- plot_enrichment_counts(
  pd_with_metadata,
  group_cols = "treatment" # makes a separate panel for each group in this metadata column
)
#> [19:51:43] INFO  Plotting enrichment counts (<phip_data>)
#>                  -> group_cols: 'treatment'
#> [19:51:43] INFO  Full-cross detected; pruning non-existing rows before plotting
#>                    - rule: keep exist == 1
#>                    - estimated reduction: ~4.3x
#> [19:51:43] INFO  building enrichment count plot
#>                  -> grouping variable: 'treatment'
#> [19:51:43] OK    plot built
#> [19:51:43] OK    building enrichment count plot - done
#>                  -> elapsed: 0.305s
#> [19:51:43] OK    Plotting enrichment counts (<phip_data>) - done
#>                  -> elapsed: 0.31s
```

``` r

p_enrichment_counts
```

![](typical-workflow_files/figure-html/plot_enriched_counts-1.png)

Now let’s calculate peptide alpha diversities in each sample. We can
also compute global (Kruskal-Wallis) and pair-wise (Wilcoxon)
comparisons and display all of this in a boxplot.

``` r

alpha_div <- compute_alpha(
  pd_with_metadata,
  group_cols = "treatment"
)
#> [19:51:44] INFO  Full-cross detected; pruning non-existing rows before alpha
#>                  calc
#>                    - rule: keep exist == 1
#>                    - estimated reduction: ~4.3x
#> [19:51:44] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'treatment'; ranks: 'peptide_id'
#> [19:51:44] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.248s

alpha_sig <- compute_alpha_significance(alpha_div)

p_alpha_richness <- plot_alpha(
  alpha_div,
  metric = "richness",
  group_col = "treatment",
  significance = alpha_sig,
  show_significance = TRUE # needs the package ggsignif
)
#> [19:51:44] INFO  plotting alpha diversity (precomputed)
#>                  -> metric: richness
#> [19:51:44] OK    plotting alpha diversity (precomputed) - done
#>                  -> elapsed: 0.118s
```

``` r

p_alpha_richness
```

![](typical-workflow_files/figure-html/display_alpha-1.png)

## Beta diversity

Please refer to the vignette [Beta Diversity Analysis in
PhIP-seq](https://polymerase3.github.io/phiper/articles/beta-diversity.html)
for detailed information.

Next, let’s calculate between-sample (beta) diversities; Jaccard
dissimilarities in this case. Based on these dissimilarities, we can
easily calculate a PCoA and a PERMANOVA.

``` r

beta_dist <- compute_distance(
  pd_with_metadata,
  distance = "Jaccard"
)
#> [19:51:44] INFO  auto-detected `value_col = "exist"` from `ps`.
#> [19:51:44] INFO  building abundance matrix from `ps` using `exist`.
#> [19:51:44] INFO  building pivot spec (sample_id x peptide_id).
#> [19:51:45] INFO  Collecting long table (sample_id, peptide_id, value).
#>                  -> compute_distance
#> [19:51:45] INFO  Pivoting to wide abundance matrix in R.
#>                  -> compute_distance
#> [19:51:45] INFO  abundance matrix has 40 samples and 2225 features after
#>                  preprocessing.
#> [19:51:45] INFO  auto normalization selected -> using none
#> [19:51:45] INFO  computing distance: jaccard
#> [19:51:46] INFO  distance matrix computation complete.

pcoa <- compute_pcoa(beta_dist)
#> [19:51:46] INFO  performing principal coordinates analysis
#> [19:51:46] INFO  extracting sample coordinates.
#> [19:51:46] INFO  summarizing eigenvalues and variance explained.
#> [19:51:46] INFO  pcoa analysis complete.

# add the treatment column to the coordinate dataframe in order to plot centroids later
pcoa$sample_coords <- pcoa$sample_coords %>%
  left_join(
    pd_with_metadata$data_long %>% select(sample_id, treatment) %>% distinct(),
    by = "sample_id",
    copy = TRUE
  )

permanova <- compute_permanova(
  beta_dist,
  ps = pd_with_metadata,
  group_col = "treatment",
  subject_col = "subject_id"
)
#> [19:51:46] INFO  preparing distance labels and metadata.
#> [19:51:46] INFO  building metadata from `ps`.
#> [19:51:46] INFO  filtering samples with missing grouping variables.
#> [19:51:46] INFO  subsetting distance matrix to complete cases.
#> [19:51:46] INFO  preparing global permanova model.
#> [19:51:46] INFO  running global permanova
#>                    - model: d_resp ~ treatment
#> [19:51:46] INFO  running pairwise permanova contrasts.

permanova_string <- paste0(
  "PERMANOVA (", permanova$n_perm, " permutations):",
  "\nF = ", round(permanova$F_stat, digits = 2),
  "; R2 = ", round(permanova$R2, digits = 2),
  "; P = ", permanova$p_value
)

p_pcoa <- plot_pcoa(
  pcoa,
  axes = c(1, 2),
  group_col = "treatment",
  ellipse_by = "group",
  show_centroids = TRUE,
  point_size = 4
) +
  labs(subtitle = permanova_string)
#> [19:51:46] INFO  Plotting PCoA: n=40 | group_col=treatment | time_col=<none> |
#>                  centroid_by=group
#>                  -> plot_pcoa
```

``` r

p_pcoa
```

![](typical-workflow_files/figure-html/display_beta-1.png)

## POP

Please refer to the vignette [Prevalence of Presence (POP) Analysis in
PhIP-seq](https://polymerase3.github.io/phiper/articles/pop-analysis.html)
for detailed information.

Comparing prevalences of peptides between two different groups is a
classic way of looking at PhIP-seq data. It’s good for exploration but
please refer to the [DELTA section](#delta) for a proper analysis of
shifts in prevalence using `phiper`’s DELTA module.

You can run POP in different (taxonomic) ranks (see peptide annotation
in the `peplib` variable). For any other ranks than `peptide_id`,
results will be aggregated by that rank. As an example we here run POP
on `peptide_id` and `species`.

``` r

pop <- compute_pop(
  x = pd_with_metadata,
  rank_cols = c("peptide_id", "species"),
  group_cols = "treatment"
) %>% tibble()
#> [19:51:46] INFO  compute_pop
#> [19:51:46] INFO  compute_pop
#>                    - ranks : peptide_id, species
#>                    - group_cols: treatment
#>                    - exist_col : exist
#>                    - pop_k_min : 1
#>                    - paired : FALSE
#> [19:51:47] INFO  ranks resolved
#>                    - available: peptide_id, species
#> [19:51:47] INFO  computing cohort sizes and validating binary group_cols
#> [19:51:48] INFO  computing presence per sample via k-of-n rule
#> [19:51:48] INFO  counting present samples per feature (pop, unpaired)
#> [19:51:48] INFO  building pairwise comparisons
#> [19:51:50] OK    materialized; computing Fisher p-values
#>                    - table: ph_pop_20260706_195148
#> [19:51:51] OK    done (compute_pop, unpaired)
#>                    - rows : 2605
#>                    - ranks : peptide_id, species
#>                    - k_min : 1
#> [19:51:51] OK    compute_pop - done
#>                  -> elapsed: 4.659s
```

Since we have POP data for two different ranks, we can make two separate
plots.

``` r

p_pop_static <- list()
for (rank_name in c("peptide_id", "species")) {

  p_pop_static[[rank_name]] <- scatter_static(
    df = pop,
    rank = rank_name,
    xlab = paste("Prevalence in", unique(pop$group1), "%"),
    ylab = paste("Prevalence in", unique(pop$group2), "%")
  ) +
    labs(title = rank_name)

}
```

You can also create interactive plots with plotly that you can later
save as HTML files.

``` r

p_pop_interactive <- list()
for (rank_name in c("peptide_id", "species")) {

  p_inter_temp <- scatter_interactive(
    df = pop,
    rank = rank_name,
    xlab = paste("Prevalence in", unique(pop$group1), "(%)"),
    ylab = paste("Prevalence in", unique(pop$group2), "(%)"),
    peplib = peplib
  )

  p_pop_interactive[[rank_name]] <- plotly::layout(
    p_inter_temp,
    autosize = TRUE,
    margin   = list(l = 70, r = 30, t = 50, b = 50),
    xaxis    = list(range = c(-2, 102), automargin = TRUE),
    yaxis    = list(range = c(-2, 102), automargin = TRUE),
    title = rank_name
  )

}
```

Let’s display one of the plots as an example:

``` r

p_pop_interactive$species
```

## DELTA

Please refer to the vignette [Delta Analysis in
PhIP-seq](https://polymerase3.github.io/phiper/articles/delta-analysis.html)
for detailed information.

The DELTA module can be computationally quite expensive, so if you have
a large dataset you may want to use the `future` package for
multithreading and run it on a high-performance computer. `future`
doesn’t work in interactive R Studio sessions, though, so we will run it
on one core here.

``` r

delta <- compute_delta(
  x = pd_with_metadata,
  rank_cols = c("peptide_id", "species"),
  group_cols = "treatment",
  B_permutations = 100
) %>%
  filter(m_eff > 5) # highly recommended: results based on fewer than 5 peptides can be unreliable
```

The most common way to visualize these results is by a forest plot to
display the features whose prevalence shifted the most.

``` r

p_delta_static <- list()
p_delta_interactive <- list()
for (rank_name in c("peptide_id", "species")) {

  p_delta_static[[rank_name]] <- forestplot(
    results_tbl = delta,
    rank_of_interest = rank_name,
    use_diverging_colors = TRUE,
    filter_significant = "p_perm",
    left_label = paste0("More in ", unique(delta$group1)),
    right_label = paste0("More in ", unique(delta$group2))
  )

  p_delta_interactive[[rank_name]] <- forestplot_interactive(
    results_tbl = delta,
    rank_of_interest = rank_name,
    use_diverging_colors = TRUE,
    filter_significant = "p_perm",
    left_label = paste0("More in ", unique(delta$group1)),
    right_label = paste0("More in ", unique(delta$group2))
  )

}
```

Let’s again look at one of these plots as an example:

``` r

p_delta_interactive$species$plot
```

## Save output

We assume you already know how to save tables and ggplot objects to
disk. The interactive plots can for example be saved to HTML with
[`htmlwidgets::saveWidget()`](https://rdrr.io/pkg/htmlwidgets/man/saveWidget.html).

`phip_data` objects can be stored as `.parquet` files using `phiperio`’s
[`export_parquet()`](https://polymerase3.github.io/phiperio/reference/export_parquet.html).
