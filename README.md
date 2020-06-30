# mg_bao

This respository contains all of the scripts needed to perform the analysis and
create the plots from Pardo & Spergel (2020) (link will be included when
available).

The [`mg_bao/`](mg_bao/) directory contains the main calculation scripts.

The [`drivers/`](drivers/) directory contains the scripts that are actually run to make all
of the plots and perform the analyses.

The [`data/`](data/) directory contains the data I used for the analysis (see the paper
for all references).

To run the scripts for yourself, create a 'results' directory and a
`data_products` directory within it. Make sure that the `RERUN_ANALYSIS` flag
within [`mg_bao.constants`](mg_bao/constants.py) is set to `True`. 
Then, run the [`plot_all.sh`](drivers/plot_all.sh) shell script.