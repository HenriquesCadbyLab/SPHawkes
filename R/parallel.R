# Parallel workers for the fitting functions.
#
# The `nprocs` argument of fit_semihawkes(), em_step(), trigger_evaluate()
# and trigger_spatial_marks() takes a whole number of workers (1 runs
# serially), "auto", or a cluster made with parallel::makeCluster(). The
# workers are socket (PSOCK) clusters, which work the same way on Windows,
# macOS and Linux. Forked workers (parallel::mclapply()) would avoid copying
# data but do not exist on Windows and are unsafe in GUIs such as RStudio on
# macOS.
#
# fit_semihawkes() starts one cluster for the whole fit and passes it down,
# so the start-up cost (starting R processes, loading packages in them) is
# paid once rather than at every call.

# Number of workers "auto" uses: one per physical core but one, leaving that
# core for the operating system and other programs, though 2 on a 2-core
# machine (the main R session mostly waits while the workers run). At most 2
# while R CMD check runs, as CRAN requires, and no more than the `mc.cores`
# option when a user has set it.
.auto_workers <- function() {
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores)) cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores)) cores <- 1L
  cores <- as.integer(cores)
  workers <- max(min(2L, cores), cores - 1L)

  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  checking <- (nzchar(limit) && !identical(toupper(limit), "FALSE")) ||
    nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_", ""))
  if (checking) workers <- min(workers, 2L)

  mc_cores <- getOption("mc.cores")
  if (is.numeric(mc_cores) && length(mc_cores) == 1L && !is.na(mc_cores) &&
      mc_cores >= 1) {
    workers <- min(workers, as.integer(mc_cores))
  }
  as.integer(workers)
}

# Smallest number of events for which "auto" uses more than one worker.
# Below it, a fit takes a few seconds and starting the workers costs about
# as much as they save.
.auto_min_events <- 1000L

# Stops unless `nprocs` is a whole number of at least 1, "auto" or a cluster.
.check_nprocs <- function(nprocs) {
  if (inherits(nprocs, "cluster") || identical(nprocs, "auto")) {
    return(invisible(nprocs))
  }
  if (!is.numeric(nprocs) || length(nprocs) != 1L || is.na(nprocs) ||
      nprocs < 1 || nprocs != round(nprocs)) {
    stop("`nprocs` must be a whole number of workers (1 runs serially), ",
         "\"auto\", or a cluster from `parallel::makeCluster()`.", call. = FALSE)
  }
  invisible(nprocs)
}

# Resolves `nprocs` to either a cluster (returned as given) or a number of
# workers. `n`, the number of events, lets "auto" stay serial for small data.
.resolve_workers <- function(nprocs, n = NULL) {
  .check_nprocs(nprocs)
  if (inherits(nprocs, "cluster")) {
    return(nprocs)
  }
  if (identical(nprocs, "auto")) {
    if (!is.null(n) && n < .auto_min_events) return(1L)
    return(.auto_workers())
  }
  cores <- parallel::detectCores()
  if (!is.na(cores) && nprocs > cores) {
    warning(sprintf(paste0("`nprocs = %d` is more than the %d cores this ",
                           "machine has, which would only slow it down; ",
                           "using %d."), as.integer(nprocs), cores, cores),
            call. = FALSE)
    nprocs <- cores
  }
  as.integer(nprocs)
}

# Starts a socket cluster of `workers` R processes and gives them this
# session's library paths, so that they load the same packages (for
# instance from an renv project library) rather than R's defaults.
#
# The sockets are opened with TCP_NODELAY ("no-delay") at both ends.
# Without it, the operating system holds back the end of each result a
# worker sends until the main session acknowledges the start, which it
# delays in turn: about 20 ms per result on Linux, which added a second to
# every trigger_evaluate() call on the Milwaukee example.
.make_cluster <- function(workers) {
  old <- options(socketOptions = "no-delay")
  on.exit(options(old), add = TRUE)
  cl <- parallel::makeCluster(
    workers,
    rscript_args = c("-e", shQuote("options(socketOptions = 'no-delay')"))
  )
  set_paths <- function(paths) invisible(.libPaths(paths))
  # A function from base only, so the workers can run it before they can
  # find this package.
  environment(set_paths) <- baseenv()
  parallel::clusterCall(cl, set_paths, .libPaths())
  cl
}

# Uses the cluster `nprocs` if it is one, or starts `nprocs` workers if that
# is more than one; returns NULL for serial work. `owned` is TRUE when the
# caller started the cluster here and must stop it.
.open_workers <- function(nprocs, n = NULL) {
  workers <- .resolve_workers(nprocs, n)
  if (inherits(workers, "cluster")) {
    return(list(cl = workers, owned = FALSE))
  }
  if (workers > 1L) {
    return(list(cl = .make_cluster(workers), owned = TRUE))
  }
  list(cl = NULL, owned = FALSE)
}

# Splits `n` items into chunks for `.run_chunks()`. The chunks depend only on
# `n`, not on the number of workers, and results are added up in chunk order,
# so a result is the same to the last bit with any number of workers. About
# 64 items per chunk, and at most 64 chunks, so that up to 16 workers get
# several chunks each to balance the load between fast and slow cores.
.fixed_chunks <- function(n, size = 64L, max_chunks = 64L) {
  parallel::splitIndices(n, min(max_chunks, max(1L, ceiling(n / size))))
}

# Runs `fun(chunk, data)` for every chunk in `chunks`, in parallel on the
# cluster `cl` if there is one, and returns the results in chunk order.
# `data` is sent to each worker once, rather than with every chunk, and
# removed again afterwards. `fun` must use only its arguments and functions
# from other packages, as it runs in the workers.
.run_chunks <- function(cl, chunks, fun, data) {
  if (is.null(cl)) {
    return(lapply(chunks, fun, data))
  }
  staged <- new.env(parent = emptyenv())
  staged$.spHawkes_chunk_work <- list(fun = fun, data = data)
  parallel::clusterExport(cl, ".spHawkes_chunk_work", envir = staged)
  clear <- function() {
    rm(list = ".spHawkes_chunk_work", envir = globalenv())
    invisible(NULL)
  }
  environment(clear) <- baseenv()
  on.exit(parallel::clusterCall(cl, clear), add = TRUE)

  run <- function(chunk) {
    work <- get(".spHawkes_chunk_work", envir = globalenv())
    work$fun(chunk, work$data)
  }
  environment(run) <- baseenv()
  parallel::clusterApplyLB(cl, chunks, run)
}
