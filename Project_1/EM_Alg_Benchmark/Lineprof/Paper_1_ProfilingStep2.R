library(lineprof)

source("~/source_EMAlgorithm_code.R")

profile <- lineprof(OverallAGAlg(dat$X, dat$Y, Bres, m_seq = 3, tuningpA = 0))

saveRDS(profile, file = "Profiling_Benchmark0.rds")

shine(profile)