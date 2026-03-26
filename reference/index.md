# Package index

## All functions

- [`add_exist()`](https://polymerase3.github.io/phiper/reference/add_exist.md)
  :

  Ensure an existence flag (all ones) on `data_long`

- [`compute_alpha_diversity()`](https://polymerase3.github.io/phiper/reference/compute_alpha_diversity.md)
  : Compute alpha diversity per sample / group across ranks

- [`compute_capscale()`](https://polymerase3.github.io/phiper/reference/compute_capscale.md)
  : Constrained Ordination (db-rda / cap) on Distance Matrix

- [`compute_delta()`](https://polymerase3.github.io/phiper/reference/compute_delta.md)
  : Global Shift in Peptide-level Prevalence via Subject-level
  Permutation

- [`compute_dispersion()`](https://polymerase3.github.io/phiper/reference/compute_dispersion.md)
  : Test Homogeneity of Dispersion (Beta Dispersion)

- [`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md)
  : Compute Pairwise Sample Distances

- [`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md)
  : Principal Components Analysis (PCoA) on a Distance Matrix

- [`compute_pcoa_feature_associations()`](https://polymerase3.github.io/phiper/reference/compute_pcoa_feature_associations.md)
  : Compute Feature Associations to PCoA Vectors

- [`compute_permanova()`](https://polymerase3.github.io/phiper/reference/compute_permanova.md)
  : PERMANOVA with Global and Post-hoc Tests on Beta Diversity

- [`compute_tsne()`](https://polymerase3.github.io/phiper/reference/compute_tsne.md)
  : Compute t-SNE Embeddings for Sample Distances

- [`deltaplot()`](https://polymerase3.github.io/phiper/reference/deltaplot.md)
  : Delta-prevalence vs Pooled Prevalence

- [`deltaplot_interactive()`](https://polymerase3.github.io/phiper/reference/deltaplot_interactive.md)
  : Interactive Delta-prevalence vs Pooled Prevalence

- [`disconnect()`](https://polymerase3.github.io/phiper/reference/disconnect.md)
  : Disconnect backend database connections

- [`.as_prev_result()`](https://polymerase3.github.io/phiper/reference/dot-as_prev_result.md)
  :

  build a `ph_prev_result` object

- [`ecdf_plot()`](https://polymerase3.github.io/phiper/reference/ecdf_plot.md)
  : ECDF of Per-peptide Prevalences for Two Groups

- [`ecdf_plot_interactive()`](https://polymerase3.github.io/phiper/reference/ecdf_plot_interactive.md)
  : ECDF of Per-peptide Prevalence for Two Groups

- [`expand_phip_data()`](https://polymerase3.github.io/phiper/reference/expand_phip_data.md)
  :

  Expand to a full `sample_id * peptide_id` grid

- [`export_parquet()`](https://polymerase3.github.io/phiper/reference/export_parquet.md)
  : Export a phip_data Table to Parquet

- [`forestplot()`](https://polymerase3.github.io/phiper/reference/forestplot.md)
  : Forest Plot of Top/Bottom Raw Stouffer T by Rank

- [`forestplot_interactive()`](https://polymerase3.github.io/phiper/reference/forestplot_interactive.md)
  : Interactive Forest Plot of Top/Bottom DELTA/Stouffer Statistics

- [`get_comparisons()`](https://polymerase3.github.io/phiper/reference/get_comparisons.md)
  : Retrieve the comparisons definition table

- [`get_counts()`](https://polymerase3.github.io/phiper/reference/get_counts.md)
  : Retrieve the main PhIP-Seq counts table

- [`get_meta()`](https://polymerase3.github.io/phiper/reference/get_meta.md)
  : Retrieve the metadata list

- [`get_peptide_library()`](https://polymerase3.github.io/phiper/reference/get_peptide_library.md)
  : Retrieve the peptide-library annotation table

- [`get_peptide_meta()`](https://polymerase3.github.io/phiper/reference/get_peptide_meta.md)
  : Retrieve the peptide metadata table into DuckDB, forcing atomic
  types

- [`merge(`*`<phip_data>`*`)`](https://polymerase3.github.io/phiper/reference/merge.phip_data.md)
  :

  Merge or join a `phip_data` object

- [`new_phip_data()`](https://polymerase3.github.io/phiper/reference/new_phip_data.md)
  :

  Construct a **phip_data** object

- [`ph_prevalence_compare()`](https://polymerase3.github.io/phiper/reference/ph_prevalence_compare.md)
  : Prevalence by group with pairwise tests (POP; per-rank FDR; BH &
  weighted BH)

- [`phip_convert()`](https://polymerase3.github.io/phiper/reference/phip_convert.md)
  :

  Convert raw PhIP-Seq output into a `phip_data` object

- [`phip_convert_legacy()`](https://polymerase3.github.io/phiper/reference/phip_convert_legacy.md)
  :

  Convert legacy Carlos-style input to a modern **phip_data** object

- [`left_join(`*`<phip_data>`*`)`](https://polymerase3.github.io/phiper/reference/phip_data_join.md)
  [`right_join(`*`<phip_data>`*`)`](https://polymerase3.github.io/phiper/reference/phip_data_join.md)
  [`inner_join(`*`<phip_data>`*`)`](https://polymerase3.github.io/phiper/reference/phip_data_join.md)
  [`full_join(`*`<phip_data>`*`)`](https://polymerase3.github.io/phiper/reference/phip_data_join.md)
  [`semi_join(`*`<phip_data>`*`)`](https://polymerase3.github.io/phiper/reference/phip_data_join.md)
  [`anti_join(`*`<phip_data>`*`)`](https://polymerase3.github.io/phiper/reference/phip_data_join.md)
  :

  dplyr joins for `phip_data`

- [`phip_example_path()`](https://polymerase3.github.io/phiper/reference/phip_example_path.md)
  : Path to example PhIP-Seq datasets shipped with phiper

- [`phip_load_example_data()`](https://polymerase3.github.io/phiper/reference/phip_load_example_data.md)
  : Load Example PhIP-Seq Dataset as \<phip_data\>

- [`phip_palette`](https://polymerase3.github.io/phiper/reference/phip_palette.md)
  : PHIP default colour palette

- [`phip_use_montserrat()`](https://polymerase3.github.io/phiper/reference/phip_use_montserrat.md)
  : Register and enable the Montserrat font for plotting (showtext)

- [`plot_alpha_diversity()`](https://polymerase3.github.io/phiper/reference/plot_alpha_diversity.md)
  : Plot alpha diversity (richness/Shannon/Simpson) from precomputed
  results

- [`plot_alpha_diversity_interactive()`](https://polymerase3.github.io/phiper/reference/plot_alpha_diversity_interactive.md)
  : Plot alpha diversity (precomputed) — interactive (plotly)

- [`plot_cap()`](https://polymerase3.github.io/phiper/reference/plot_cap.md)
  : Plot CAP/db-RDA Results (Constrained Ordination)

- [`plot_dispersion()`](https://polymerase3.github.io/phiper/reference/plot_dispersion.md)
  : Plot beta dispersion (distance to centroid)

- [`plot_enrichment_counts()`](https://polymerase3.github.io/phiper/reference/plot_enrichment_counts.md)
  : Plot enrichment counts per group (and optional interaction)

- [`plot_pcoa()`](https://polymerase3.github.io/phiper/reference/plot_pcoa.md)
  : Plot Principal Coordinates Analysis (PCoA)

- [`plot_scree()`](https://polymerase3.github.io/phiper/reference/plot_scree.md)
  : Scree Plot for PCoA Eigenvalues

- [`plot_tsne()`](https://polymerase3.github.io/phiper/reference/plot_tsne.md)
  : Plot t-SNE embeddings

- [`prev_filter_pairs()`](https://polymerase3.github.io/phiper/reference/prev_filter_pairs.md)
  : Filter pairwise results by groups/ranks/features with optional
  q-value gates

- [`scale_colour_phip()`](https://polymerase3.github.io/phiper/reference/scale_colour_phip.md)
  [`scale_color_phip()`](https://polymerase3.github.io/phiper/reference/scale_colour_phip.md)
  [`scale_fill_phip()`](https://polymerase3.github.io/phiper/reference/scale_colour_phip.md)
  : Discrete colour & fill scales using the PHIP palette

- [`scatter_interactive()`](https://polymerase3.github.io/phiper/reference/scatter_interactive.md)
  :

  Interactive prevalence scatter for `ph_prev_result`

- [`scatter_static()`](https://polymerase3.github.io/phiper/reference/scatter_static.md)
  :

  Static scatterplot of percent1 vs percent2 from
  [`ph_prevalence_compare()`](https://polymerase3.github.io/phiper/reference/ph_prevalence_compare.md)

- [`summary(`*`<ph_prev_result>`*`)`](https://polymerase3.github.io/phiper/reference/summary.ph_prev_result.md)
  :

  Summarize FDR accounting and design for `ph_prev_result`

- [`theme_phip()`](https://polymerase3.github.io/phiper/reference/theme_phip.md)
  :

  Theme `theme_phip`

- [`volcano_interactive()`](https://polymerase3.github.io/phiper/reference/volcano_interactive.md)
  : Interactive volcano plot (log2 ratio vs -log10 p)

- [`volcano_static()`](https://polymerase3.github.io/phiper/reference/volcano_static.md)
  : Static volcano plot (log2 ratio vs -log10 p)

- [`write_result()`](https://polymerase3.github.io/phiper/reference/write_result.md)
  : Generic writer for phiper result objects

- [`write_result(`*`<ph_prev_result>`*`)`](https://polymerase3.github.io/phiper/reference/write_result.ph_prev_result.md)
  : Write method for ph_prev_result (ph_prevalence_compare output)
