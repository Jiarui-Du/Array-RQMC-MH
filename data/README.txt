Data are organized by dataset:
  ripley/author/{myX.txt,Y.txt}       exact author matrices used by the model
  ripley/source/{synth.tr.csv,source_match.json}  source data and row map
  pima/author/{myX.txt,Y.txt}         exact author matrices used by the model
  pima/source/{Pima.tr.csv,Pima.te.csv,source_match.json}  source data and row map

The author text matrices are read by load_bayesian_model.m.
Source: https://gitlab.com/aggelisalexopoulos/variance-reduction
Pinned commit: 07c5bc69c664904881f4845d4628d2ce0b8e1981
preprocess_bayesian_data.m reconstructs and verifies them from source/ CSVs.
See source_manifest.json for original URLs and hashes. Relocation changes only
local paths; source URLs, file contents, rounding and row maps are unchanged.
Ripley: 250 observations, 3 coefficients. Pima: 532 observations, 8 coefficients.
Pima.tr and Pima.te are combined to reproduce the author's 532-row experiment;
they are not separate training and held-out evaluation sets in this experiment.
The loader moves the intercept to column 1 and preserves author text rounding.
The Bayesian experiments use the original improper flat prior and start at MLE.
