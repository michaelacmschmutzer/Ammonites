#!/bin/bash

# Little bash script to run entire analysis when needed

# Prelude. Compile and clean dataset and make some preliminary plots
Rscript merge_datasets.R
Rscript make_diagnostic_plots.R
Rscript make_density_maps.R

# Toccata. Run the main body of GIS analyses
Rscript rotate_positions.R
Rscript draw_convex_hulls.R 'ammonoids'
Rscript draw_convex_hulls.R 'nautilids'
Rscript subsample_distributions.R 'ammonoids'
Rscript subsample_distributions.R 'nautilids'
Rscript impose_palaeocoasts.R 'ammonoids'
Rscript impose_palaeocoasts.R 'nautilids'

# Fugue. Other analyses, with stats
Rscript plot_hatching_size.R &>/dev/null
Rscript plot_body_size.R &>/dev/null
Rscript analyse_abundance.R &>/dev/null
Rscript analyse_subsampling.R 'ammonoids' &>/dev/null
Rscript analyse_subsampling.R 'nautilids' &>/dev/null
RScript compare_hypotheses.R &>/dev/null

# Finale. Figure plotting
