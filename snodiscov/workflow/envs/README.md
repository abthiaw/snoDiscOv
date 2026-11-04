Conda environments (channels: conda-forge, bioconda, nodefaults; strict priority
recommended: `conda config --set channel_priority strict`).

Versions come from the environments used for the article. Packages that were not
pinned there (multiqc, pandas/numpy in python.yaml, r-tidyverse, python in
coco.yaml) carry a version chosen for compatibility: confirm them by exporting
the environments actually used (`snodiscov export-envs`, see the main README),
which also writes `<env>.linux-64.pin.txt` files. Snakemake uses these pin files
automatically when they sit next to the .yaml, giving bit-for-bit identical
environments.

The `anaconda`/`defaults` channels used originally were removed: they are subject
to Anaconda's terms of service and mixing them with conda-forge hurts reproducibility.
