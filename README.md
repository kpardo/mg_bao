# mg_bao

This respository contains all of the scripts needed to perform the analysis and
create the plots from [Pardo & Spergel (2020)](https://ui.adsabs.harvard.edu/abs/2020arXiv200700555P/abstract).

The [`mg_bao/`](mg_bao/) directory contains the main calculation scripts.

The [`drivers/`](drivers/) directory contains the scripts that are actually run to make all
of the plots and perform the analyses.

The [`data/`](data/) directory contains the data I used for the analysis (see the paper
for all references).

To run the scripts for yourself, make sure that the `RERUN_ANALYSIS` flag
within [`mg_bao.constants`](mg_bao/constants.py) is set to `True`. 
Then, run the [`plot_all.sh`](drivers/plot_all.sh) shell script.
