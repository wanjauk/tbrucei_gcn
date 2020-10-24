# ensure results are reproducible
set.seed(1)

# other settings
options(digits = 4)
options(stringsAsFactors = FALSE)

# set number of threads to allow in WGCNA analysis
allowWGCNAThreads(nThreads=20)