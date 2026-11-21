#!/usr/bin/env Rscript
# Project per-step walltime, CPU, memory, SUs and disk for a new euk-bot run from the raw read sizes of its samples.
#
# Calibration comes from finished runs: Nextflow trace files (realtime, %cpu, peak_rss, wchar) joined to sacct
# (allocated CPUs, requested memory, elapsed, time limit), and the run's qc/preliminary_qc.tsv (raw_pairs per sample).
# New samples are sized from their FASTQ file sizes (bytes per read pair taken from the calibration samples),
# or from a raw_pairs column in the new samplesheet if there is one.
#
# Scaling model (see SCALING below): metric_new = metric_obs * (size_new / size_obs)^exponent, where size is the
# sample's raw pairs for per-sample steps and the run's total raw pairs for pooled steps (clustering, annotation).
# Each per-sample step gives a "typical" projection (median over calibration samples) and a "high" one (the
# calibration sample that needs the most per pair). Only steps seen in a calibration run can be projected;
# re-run this after a run that exercises new branches.
#
# Usage (on Grace: module load GCC/14.3.0 R/4.5.2; light enough for a login node, it only reads tables and calls sacct):
#   Rscript tools/project_resources.R \
#       --calib /scratch/user/skhu/euk-bot/test3_v2 \
#       --calib-samplesheet samplesheet_test3.csv \
#       --samplesheet samplesheet_full.csv \
#       --out projection
#   --calib may be given several times (comma-separated) for several finished runs; the newest task per step wins.
#   Writes <out>.steps.tsv (one row per step) and <out>.tasks.tsv (one row per projected task), prints a summary.
#
# Options:
#   --mem-per-core-gb 7.5   Grace charges max(cores, ceil(mem / this)) core-hours as SUs; check HPRC's current policy
#   --node-mem-gb 360       usable memory on a standard node (flag tasks that cannot fit)
#   --max-walltime-h 168    longest partition walltime used (long = 7 d)
#   --sacct-file FILE       saved `sacct -X -P -n --format=JobID,AllocCPUS,ReqMem,ElapsedRaw,TimelimitRaw,State`
#                           output, for when sacct no longer has old jobs

suppressWarnings(suppressMessages(library(utils)))

# ---------------------------------------------------------------- arguments
parse_args <- function(argv) {
  opts <- list(calib = character(), calib_samplesheet = character(), samplesheet = NULL, out = "projection",
               mem_per_core_gb = 7.5, node_mem_gb = 360, max_walltime_h = 168, sacct_file = NULL)
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i]); key <- gsub("-", "_", key)
    if (i == length(argv)) stop("missing value for --", key)
    val <- argv[i + 1]
    if (key %in% c("calib", "calib_samplesheet")) {
      opts[[key]] <- c(opts[[key]], strsplit(val, ",")[[1]])
    } else if (key %in% names(opts)) {
      opts[[key]] <- if (is.numeric(opts[[key]])) as.numeric(val) else val
    } else stop("unknown option --", key)
    i <- i + 2
  }
  if (!length(opts$calib)) stop("give --calib <finished run outdir>")
  if (!length(opts$calib_samplesheet)) stop("give --calib-samplesheet (the samplesheet of the calibration run)")
  if (is.null(opts$samplesheet)) stop("give --samplesheet (the run to project)")
  opts
}

# ---------------------------------------------------------------- parsers for Nextflow trace fields
parse_duration_h <- function(x) {
  vapply(x, function(s) {
    if (is.na(s) || s %in% c("-", "")) return(NA_real_)
    parts <- regmatches(s, gregexpr("[0-9.]+(ms|d|h|m|s)", s))[[1]]
    if (!length(parts)) return(NA_real_)
    mult <- c(d = 24, h = 1, m = 1 / 60, s = 1 / 3600, ms = 1 / 3.6e6)
    sum(vapply(parts, function(p) {
      u <- sub("^[0-9.]+", "", p); as.numeric(sub("(ms|d|h|m|s)$", "", p)) * mult[[u]]
    }, numeric(1)))
  }, numeric(1), USE.NAMES = FALSE)
}

parse_size_gb <- function(x) {
  mult <- c(B = 1 / 1024^3, KB = 1 / 1024^2, MB = 1 / 1024, GB = 1, TB = 1024, PB = 1024^2)
  vapply(x, function(s) {
    if (is.na(s) || s %in% c("-", "")) return(NA_real_)
    m <- regmatches(s, regexec("^([0-9.]+) *([KMGTP]?B)$", trimws(s)))[[1]]
    if (length(m) != 3) return(NA_real_)
    as.numeric(m[2]) * mult[[m[3]]]
  }, numeric(1), USE.NAMES = FALSE)
}

parse_reqmem_gb <- function(x) {   # sacct ReqMem: 96G, 16000M, 4Gn / 4Gc (older Slurm)
  vapply(x, function(s) {
    m <- regmatches(s, regexec("^([0-9.]+)([KMGT])", s))[[1]]
    if (length(m) != 3) return(NA_real_)
    as.numeric(m[2]) * c(K = 1 / 1024^2, M = 1 / 1024, G = 1, T = 1024)[[m[3]]]
  }, numeric(1), USE.NAMES = FALSE)
}

# ---------------------------------------------------------------- scaling assumptions per step
# time_exp / mem_exp / disk_exp: exponent on the size ratio. 0 = does not grow with input (e.g. memory set by a
# reference DB index), 1 = linear. Assemblers and k-mer tables are given sublinear memory growth. These are
# assumptions, not fits: with 3 calibration samples there is too little spread to fit exponents. Edit as data come in.
SCALING <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
step                    driver   time_exp mem_exp disk_exp
CAT_FASTQ               sample   1        0       1
FASTQC_RAW              sample   1        0       1
FASTQC_CLEAN            sample   1        0       1
FASTQC_NORM             sample   1        0       1
FASTP                   sample   1        0       1
BBDUK                   sample   1        0       1
RIBODETECTOR            sample   1        0       1
BBNORM                  sample   1        0.5     1
STRANDEDNESS            sample   0        0       0
SOURMASH_SKETCH         sample   1        0       1
NONPAREIL               sample   1        0       1
KAIJU                   sample   1        0       1
DIAMOND_BLASTX          sample   1        0       1
RNASPADES               sample   1        0.5     1
MEGAHIT                 sample   1        0.7     1
SEQKIT_STATS            sample   1        0       1
ASSEMBLY_MAPRATE        sample   1        0.5     1
BUSCO                   sample   0.5      0       0.5
SALMON_QUANT            sample   1        0       1
MULTIQC_RAW             nsample  1        1       1
MULTIQC_CLEAN           nsample  1        1       1
QC_SUMMARY              nsample  1        1       1
SOURMASH_COMPARE        nsample2 1        1       1
READ_PROFILE_SUMMARY    nsample  1        1       1
SUMMARIZE_HITS          total    1        0.5     1
ASSEMBLY_QC_SUMMARY     nsample  1        1       1
CAT_GROUP_READS         total    1        0       1
MMSEQS_LINCLUST         total    1        1       1
TRANSDECODER            total    1        0.5     1
METAEUK                 total    1        0.5     1
DIAMOND_BLASTP          total    1        0       1
EGGNOG_MAPPER           total    1        0       1
EUKULELE                total    1        0.5     1
SALMON_INDEX            total    1        1       1
MERGE_SALMON            nsample  1        1       1
CONTIG_TABLE            total    1        1       1
")
DEFAULT_SCALING <- list(driver = "total", time_exp = 1, mem_exp = 1, disk_exp = 1)  # conservative for unknown steps

# Main dependency chain, for a rough wall-clock estimate (queue wait not included). Steps on one line run in parallel.
CRITICAL_PATH <- list(
  c("FASTP"), c("BBDUK"), c("RIBODETECTOR"), c("BBNORM"), c("RNASPADES", "MEGAHIT"), c("MMSEQS_LINCLUST"),
  c("TRANSDECODER", "METAEUK"), c("DIAMOND_BLASTP"), c("EUKULELE", "EGGNOG_MAPPER"), c("CONTIG_TABLE"))

# ---------------------------------------------------------------- calibration data
read_tsv <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE, quote = "", na.strings = c("", "NA"))

read_traces <- function(outdirs) {
  files <- unlist(lapply(outdirs, function(d) Sys.glob(file.path(d, "pipeline_info", "trace_*.txt"))))
  if (!length(files)) stop("no pipeline_info/trace_*.txt under: ", paste(outdirs, collapse = ", "))
  tr <- do.call(rbind, lapply(files, function(f) {
    t <- read_tsv(f); t$trace_file <- basename(f); t$trace_order <- file.mtime(f); t
  }))
  tr <- tr[tr$status == "COMPLETED", ]
  # name like "ASSEMBLY:MEGAHIT (137)" -> step MEGAHIT, tag 137
  tr$step <- sub(" .*$", "", sub("^.*:", "", tr$name))
  tr$tag  <- ifelse(grepl("\\(", tr$name), sub("^.*\\((.*)\\)$", "\\1", tr$name), "")
  # newest completed task per step + tag (older runs may have used other code or settings)
  tr <- tr[order(tr$trace_order, decreasing = TRUE), ]
  tr <- tr[!duplicated(paste(tr$step, tr$tag)), ]
  tr$realtime_h  <- parse_duration_h(tr$realtime)
  tr$peak_rss_gb <- parse_size_gb(tr$peak_rss)
  tr$wchar_gb    <- parse_size_gb(tr$wchar)
  tr$cpu_used    <- as.numeric(sub("%", "", tr$`%cpu`)) / 100
  tr
}

read_sacct <- function(job_ids, sacct_file) {
  cols <- c("JobID", "AllocCPUS", "ReqMem", "ElapsedRaw", "TimelimitRaw", "State")
  if (!is.null(sacct_file)) {
    lines <- readLines(sacct_file)
  } else {
    if (!nzchar(Sys.which("sacct"))) { warning("sacct not found; SUs fall back to trace-only estimates"); return(NULL) }
    ids <- unique(job_ids[grepl("^[0-9]+$", job_ids)])
    lines <- unlist(lapply(split(ids, ceiling(seq_along(ids) / 200)), function(chunk)
      system2("sacct", c("-X", "-P", "-n", "-j", paste(chunk, collapse = ","),
                         paste0("--format=", paste(cols, collapse = ","))), stdout = TRUE)))
  }
  if (!length(lines)) return(NULL)
  s <- as.data.frame(do.call(rbind, strsplit(lines, "|", fixed = TRUE)), stringsAsFactors = FALSE)
  names(s) <- cols[seq_len(ncol(s))]
  data.frame(native_id = s$JobID, alloc_cpus = as.numeric(s$AllocCPUS), req_mem_gb = parse_reqmem_gb(s$ReqMem),
             elapsed_h = as.numeric(s$ElapsedRaw) / 3600, timelimit_h = as.numeric(s$TimelimitRaw) / 60,
             stringsAsFactors = FALSE)
}

read_samplesheet <- function(f) {
  ss <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("sample", "fastq_1", "fastq_2") %in% names(ss))) stop(f, ": needs columns sample,fastq_1,fastq_2")
  ss$bytes <- file.size(ss$fastq_1) + file.size(ss$fastq_2)
  if (anyNA(ss$bytes)) stop(f, ": cannot stat ", paste(ss$sample[is.na(ss$bytes)], collapse = ", "))
  agg <- aggregate(bytes ~ sample, ss, sum)   # lanes of one sample are concatenated by CAT_FASTQ
  if (!"group" %in% names(ss)) ss$group <- ""
  agg$group <- ss$group[match(agg$sample, ss$sample)]
  if ("raw_pairs" %in% names(ss)) agg$raw_pairs <- aggregate(raw_pairs ~ sample, ss, sum)$raw_pairs[match(agg$sample, unique(ss$sample))]
  agg$group[is.na(agg$group)] <- ""
  agg
}

su <- function(cores, mem_gb, hours, mem_per_core) pmax(cores, ceiling(mem_gb / mem_per_core)) * hours

fmt <- function(x, d = 1) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d, big.mark = ","))

# ---------------------------------------------------------------- main
main <- function() {
  o <- parse_args(commandArgs(trailingOnly = TRUE))

  qc <- do.call(rbind, lapply(o$calib, function(d) read_tsv(file.path(d, "qc", "preliminary_qc.tsv"))[, c("sample", "raw_pairs")]))
  qc$sample <- as.character(qc$sample); qc <- qc[!duplicated(qc$sample), ]

  calib_ss <- do.call(rbind, lapply(o$calib_samplesheet, read_samplesheet))
  calib_ss$raw_pairs <- qc$raw_pairs[match(as.character(calib_ss$sample), qc$sample)]
  calib_ss <- calib_ss[!is.na(calib_ss$raw_pairs), ]
  bytes_per_pair <- sum(calib_ss$bytes) / sum(calib_ss$raw_pairs)

  new <- read_samplesheet(o$samplesheet)
  if (is.null(new$raw_pairs)) {
    new$raw_pairs <- new$bytes / bytes_per_pair
    new$pairs_source <- "estimated from fastq.gz size"
  } else new$pairs_source <- "samplesheet raw_pairs"

  tr <- read_traces(o$calib)
  sa <- read_sacct(tr$native_id, o$sacct_file)
  if (!is.null(sa)) tr <- merge(tr, sa, by = "native_id", all.x = TRUE)
  for (col in c("alloc_cpus", "req_mem_gb", "elapsed_h", "timelimit_h")) if (is.null(tr[[col]])) tr[[col]] <- NA_real_
  tr$elapsed_h[is.na(tr$elapsed_h)] <- tr$realtime_h[is.na(tr$elapsed_h)]
  tr$alloc_cpus[is.na(tr$alloc_cpus)] <- ceiling(tr$cpu_used[is.na(tr$alloc_cpus)])

  obs_samples <- qc$sample
  obs_total <- sum(qc$raw_pairs)
  obs_n <- nrow(qc)
  new_total <- sum(new$raw_pairs)
  new_n <- nrow(new)
  tr$sample_pairs <- qc$raw_pairs[match(tr$tag, obs_samples)]

  sc_of <- function(step) {
    r <- SCALING[SCALING$step == step, ]
    if (nrow(r)) as.list(r[1, ]) else c(list(step = step), DEFAULT_SCALING)
  }

  tasks <- list(); steps <- list()
  for (key in unique(paste(tr$step, ifelse(tr$tag %in% obs_samples, "", tr$tag), sep = "\t"))) {
    st <- strsplit(key, "\t")[[1]]; step <- st[1]; label <- if (length(st) > 1) st[2] else ""
    obs <- tr[tr$step == step & (if (nzchar(label)) tr$tag == label else tr$tag %in% obs_samples | tr$tag == ""), ]
    sc <- sc_of(step)
    per_sample <- sc$driver == "sample" && all(obs$tag %in% obs_samples)
    if (sc$driver == "sample" && !per_sample) sc$driver <- "total"   # a per-sample step run on pooled input

    if (per_sample) {
      # per-pair intensity of each calibration sample; median = typical, max = high
      ref <- obs$sample_pairs
      pick <- function(v, exp) { r <- v / ref^exp; c(typical = median(r), high = max(r)) }
      t_k <- pick(obs$elapsed_h, sc$time_exp); m_k <- pick(obs$peak_rss_gb, sc$mem_exp); d_k <- pick(obs$wchar_gb, sc$disk_exp)
      targets <- data.frame(task = new$sample, size = new$raw_pairs, stringsAsFactors = FALSE)
      # co-assemblies: MEGAHIT also runs once per group on the group's pooled reads
      if (step == "MEGAHIT" && any(nzchar(new$group))) {
        g <- aggregate(raw_pairs ~ group, new[nzchar(new$group), ], sum)
        targets <- rbind(targets, data.frame(task = paste0("coassembly:", g$group), size = g$raw_pairs))
      }
      proj <- data.frame(step = step, label = label, task = targets$task, raw_pairs = targets$size,
                         hours = t_k["typical"] * targets$size^sc$time_exp,
                         hours_high = t_k["high"] * targets$size^sc$time_exp,
                         peak_mem_gb = m_k["typical"] * targets$size^sc$mem_exp,
                         peak_mem_gb_high = m_k["high"] * targets$size^sc$mem_exp,
                         disk_gb = d_k["typical"] * targets$size^sc$disk_exp, stringsAsFactors = FALSE)
      size_range <- range(ref)
    } else {
      ratio <- switch(sc$driver, nsample = new_n / obs_n, nsample2 = (new_n / obs_n)^2, new_total / obs_total)
      o1 <- obs[1, ]
      proj <- data.frame(step = step, label = label, task = if (nzchar(label)) label else "pooled",
                         raw_pairs = new_total,
                         hours = o1$elapsed_h * ratio^sc$time_exp, hours_high = o1$elapsed_h * ratio^sc$time_exp,
                         peak_mem_gb = o1$peak_rss_gb * ratio^sc$mem_exp, peak_mem_gb_high = o1$peak_rss_gb * ratio^sc$mem_exp,
                         disk_gb = o1$wchar_gb * ratio^sc$disk_exp, stringsAsFactors = FALSE)
      size_range <- rep(obs_total, 2)
    }
    cpus <- max(obs$alloc_cpus, na.rm = TRUE)
    req_mem <- suppressWarnings(max(obs$req_mem_gb, na.rm = TRUE)); if (!is.finite(req_mem)) req_mem <- NA
    tlim <- suppressWarnings(max(obs$timelimit_h, na.rm = TRUE)); if (!is.finite(tlim)) tlim <- NA
    cpu_eff <- median(obs$cpu_used / obs$alloc_cpus, na.rm = TRUE)
    proj$cpus <- cpus
    proj$req_mem_gb <- req_mem
    proj$su <- su(cpus, ifelse(is.na(req_mem), proj$peak_mem_gb, req_mem), proj$hours, o$mem_per_core_gb)
    proj$su_high <- su(cpus, ifelse(is.na(req_mem), proj$peak_mem_gb_high, pmax(req_mem, proj$peak_mem_gb_high)),
                       proj$hours_high, o$mem_per_core_gb)
    flag <- character(nrow(proj))
    add <- function(cond, msg) flag[cond] <<- ifelse(nzchar(flag[cond]), paste(flag[cond], msg, sep = "; "), msg)
    add(!is.na(tlim) & proj$hours_high > tlim, sprintf("may pass %.0f h limit (retry)", tlim))
    add(proj$hours_high > o$max_walltime_h, "longer than max partition walltime")
    add(!is.na(req_mem) & proj$peak_mem_gb_high > req_mem, sprintf("peak may pass %.0f GB request (retry)", req_mem))
    add(proj$peak_mem_gb_high > o$node_mem_gb, "needs bigmem node")
    add(proj$raw_pairs > 1.5 * max(size_range) | proj$raw_pairs < 0.5 * min(size_range), "extrapolated")
    proj$flags <- flag
    tasks[[key]] <- proj

    steps[[key]] <- data.frame(
      step = step, label = label, driver = sc$driver, n_obs = nrow(obs), n_tasks = nrow(proj),
      cpus = cpus, req_mem_gb = req_mem, timelimit_h = tlim, cpu_efficiency = cpu_eff,
      max_hours = max(proj$hours), max_hours_high = max(proj$hours_high),
      max_peak_mem_gb = max(proj$peak_mem_gb_high), total_su = sum(proj$su), total_su_high = sum(proj$su_high),
      disk_gb = sum(proj$disk_gb), time_exp = sc$time_exp, mem_exp = sc$mem_exp,
      flags = paste(unique(unlist(strsplit(flag[nzchar(flag)], "; "))), collapse = "; "), stringsAsFactors = FALSE)
  }
  tasks <- do.call(rbind, tasks); steps <- do.call(rbind, steps)
  rownames(tasks) <- NULL; rownames(steps) <- NULL
  steps <- steps[order(-steps$total_su_high), ]

  # rough wall clock: sum of the slowest task per stage on the main chain, plus 2-core head job for that long
  stage_h <- function(col) sum(vapply(CRITICAL_PATH, function(s) {
    v <- tasks[[col]][tasks$step %in% s]; if (length(v)) max(v) else 0 }, numeric(1)))
  wall <- c(typical = stage_h("hours"), high = stage_h("hours_high"))
  head_su <- 2 * wall

  write.table(steps, paste0(o$out, ".steps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tasks, paste0(o$out, ".tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- console summary
  options(width = 200)
  cat(sprintf("Calibration: %d samples, %s M raw pairs; %.1f bytes of fastq.gz per pair\n", obs_n,
              fmt(obs_total / 1e6, 1), bytes_per_pair))
  cat(sprintf("Projection : %d samples, %s M raw pairs (%s); %d co-assembly groups\n", new_n, fmt(new_total / 1e6, 1),
              unique(new$pairs_source), length(unique(new$group[nzchar(new$group)]))))
  cat(sprintf("  per sample: %s - %s M pairs\n\n", fmt(min(new$raw_pairs) / 1e6), fmt(max(new$raw_pairs) / 1e6)))
  show <- data.frame(step = ifelse(nzchar(steps$label), paste0(steps$step, " (", steps$label, ")"), steps$step),
                     tasks = steps$n_tasks, cpus = steps$cpus, req_GB = fmt(steps$req_mem_gb, 0),
                     peak_GB = fmt(steps$max_peak_mem_gb, 0), longest_h = fmt(steps$max_hours_high),
                     SUs = fmt(steps$total_su, 0), SUs_high = fmt(steps$total_su_high, 0), stringsAsFactors = FALSE)
  print(show, row.names = FALSE, right = TRUE)
  flagged <- steps[nzchar(steps$flags), ]
  if (nrow(flagged)) {
    cat("\nFlags (high projection; 'extrapolated' = new size is outside 0.5-1.5x the calibration sizes):\n")
    for (i in seq_len(nrow(flagged))) cat(sprintf("  %-32s %s\n",
      if (nzchar(flagged$label[i])) paste0(flagged$step[i], " (", flagged$label[i], ")") else flagged$step[i], flagged$flags[i]))
  }
  cat(sprintf("\nTotal SUs  : %s typical, %s high (+ head job %s-%s)\n", fmt(sum(steps$total_su), 0),
              fmt(sum(steps$total_su_high), 0), fmt(head_su[1], 0), fmt(head_su[2], 0)))
  cat(sprintf("Wall clock : ~%s-%s h on the main chain, queue wait not included\n", fmt(wall[1]), fmt(wall[2])))
  cat(sprintf("I/O        : ~%s GB written by tasks (wchar: includes temporary files; not the space needed at once)\n",
              fmt(sum(steps$disk_gb), 0)))
  cat(sprintf("SU rule    : max(cpus, ceil(GB / %.1f)) x hours, using each step's current cpus and memory request\n",
              o$mem_per_core_gb))
  # processes defined in modules/local but never seen in calibration (aliases such as FASTQC_RAW count as FASTQC)
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  mod_files <- Sys.glob(file.path(dirname(dirname(normalizePath(script))), "modules", "local", "*.nf"))
  procs <- unique(unlist(lapply(mod_files, function(f)
    sub("^process +([A-Za-z0-9_]+).*$", "\\1", grep("^process +[A-Za-z0-9_]+", readLines(f), value = TRUE)))))
  seen <- unique(tr$step)
  unseen <- procs[!vapply(procs, function(p) any(seen == p | startsWith(seen, paste0(p, "_"))), logical(1))]
  if (length(unseen)) cat("Not in calibration, so not projected:", paste(sort(unseen), collapse = ", "), "\n")
  cat(sprintf("Wrote %s.steps.tsv and %s.tasks.tsv\n", o$out, o$out))
}

main()
