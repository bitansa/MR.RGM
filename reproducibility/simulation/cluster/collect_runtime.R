###############################################################################
## collect_runtime.R  --  combines the runtime blocks of the cluster run into the files of RunTimeAnalysis.R
##   Runtime Analysis/runtime_vs_p_n30000.csv           median seconds per method and p (NA: not fitted at that p)
##   Runtime Analysis/runtime_vs_p_n30000_all_runs.csv  every timed run
##   Runtime Analysis/runtime_versions.txt              package versions, nodes and CPU
## The figure is then drawn with  RUNTIME_PLOT_ONLY <- TRUE; source("Runtime Analysis/RunTimeAnalysis.R")
###############################################################################
source("sim_functions.R")
fs <- list.files(file.path("Runtime Analysis", "cluster_parts"), "^rt_.*\\.rds$", full.names = TRUE)
if (!length(fs)) stop("no runtime blocks in Runtime Analysis/cluster_parts")
parts <- lapply(fs, readRDS)
if (!all(vapply(parts, function(x) identical(attr(x, "sampler"), SAMPLER_TAG), logical(1)))) stop("runtime block without the current sampler tag")
raw <- do.call(rbind, parts)
key <- paste(raw$Method, raw$p)
rt  <- do.call(rbind, lapply(split(raw, key), function(d) data.frame(Method = d$Method[1], p = d$p[1],
         median_sec = if (any(d$failed)) NA_real_ else median(d$sec), runs = if (any(d$failed)) 0L else sum(!is.na(d$sec)))))
rt  <- rt[order(rt$p, match(rt$Method, names(RUNTIME_EXPRS))), ]; rownames(rt) <- NULL
print(rt, row.names = FALSE)
if (any(rt$runs > 0 & rt$runs != 100)) warning("some methods have fewer than 100 timed runs")
write.csv(rt[, c("Method", "p", "median_sec")], file.path("Runtime Analysis", "runtime_vs_p_n30000.csv"), row.names = FALSE)
write.csv(raw[!raw$failed, c("Method", "p", "sec")], file.path("Runtime Analysis", "runtime_vs_p_n30000_all_runs.csv"), row.names = FALSE)
pv <- attr(parts[[1]], "package_versions")
writeLines(c(paste("Cluster run; collected on", format(Sys.time())),
             paste("CPU:", paste(unique(vapply(parts, function(x) as.character(attr(x, "cpu")), "")), collapse = "; ")),
             paste("Nodes:", paste(sort(unique(vapply(parts, function(x) attr(x, "node"), ""))), collapse = ", ")),
             "Each run on one core; the 100 runs of a method and p in blocks of 2 (MrDAG, GrIVET) or 10 runs, one block per task, each after an untimed warm-up fit.",
             paste(names(pv), pv, sep = " ")), file.path("Runtime Analysis", "runtime_versions.txt"))
