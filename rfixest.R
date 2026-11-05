# rfixest 1.2.0: R side of the Stata command rfixest. Called by rfixest.ado; not meant to be run by hand.
t_start <- proc.time()[["elapsed"]]
argv <- commandArgs(trailingOnly = TRUE)
dir <- argv[1]
mode <- if (length(argv) > 1) argv[2] else "fit"

logcon <- file(file.path(dir, "log.txt"), open = "wt")
sink(logcon)
sink(logcon, type = "message")
options(warn = 1)

fail <- function(msg) {
  writeLines(msg, file.path(dir, "err.txt"))
  cat(msg, "\n")
  sink(type = "message")
  sink()
  close(logcon)
  quit(save = "no", status = 1)
}

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

if (mode == "setup") {
  if (!requireNamespace("fixest", quietly = TRUE)) {
    lib <- .libPaths()[1]
    if (file.access(lib, 2) != 0) {
      lib <- user_lib
      dir.create(lib, recursive = TRUE, showWarnings = FALSE)
      .libPaths(c(lib, .libPaths()))
    }
    res <- try(install.packages("fixest", lib = lib, repos = "https://cloud.r-project.org"), silent = TRUE)
    if (!requireNamespace("fixest", quietly = TRUE)) fail("Could not install the R package fixest. See the R output above.")
  }
  writeLines(as.character(utils::packageVersion("fixest")), file.path(dir, "ok.txt"))
  sink(type = "message")
  sink()
  close(logcon)
  quit(save = "no", status = 0)
}

if (!requireNamespace("fixest", quietly = TRUE)) fail("RFIXEST_NO_FIXEST")

res <- tryCatch({
  meta <- strsplit(readLines(file.path(dir, "meta.txt")), " ", fixed = TRUE)
  names(meta) <- vapply(meta, `[`, "", 1L)
  meta <- lapply(meta, `[`, -1L)
  one <- function(k) meta[[k]][1L]

  suppressPackageStartupMessages(library(fixest))
  setFixest_notes(FALSE)
  t_lib <- proc.time()[["elapsed"]]

  d <- foreign::read.dta(file.path(dir, "d.dta"), convert.dates = FALSE, convert.factors = FALSE, missing.type = FALSE)
  t_read <- proc.time()[["elapsed"]]

  Y   <- meta$y
  X   <- meta$x
  FE  <- meta$fe
  CL  <- meta$cl
  W   <- meta$w
  vce <- one("vce")
  fam <- one("family")
  nth <- as.integer(one("threads"))
  tol <- as.numeric(one("tol"))
  keep_single <- one("keepsingletons") == "1"

  rhs <- if (length(X)) paste(X, collapse = "+") else "1"
  make_fml <- function(ys) {
    lhs <- if (length(ys) > 1L) paste0("c(", paste(ys, collapse = ","), ")") else ys
    f <- paste(lhs, "~", rhs)
    if (length(FE)) f <- paste(f, "|", paste(FE, collapse = "+"))
    as.formula(f, env = globalenv())
  }
  same_na <- length(Y) == 1L || all(vapply(Y[-1L], function(v) identical(is.na(d[[v]]), is.na(d[[Y[1L]]])), TRUE))
  fml <- make_fml(Y)

  vc <- switch(vce,
    cluster = as.formula(paste("~", paste(CL, collapse = "+")), env = globalenv()),
    robust  = "hetero",
    "iid")

  args <- list(fml = fml, data = d, vcov = vc, fixef.tol = tol)
  if (nth > 0) args$nthreads <- nth
  if (length(W)) args$weights <- as.formula(paste("~", W), env = globalenv())
  if (one("ppmlse") == "1") args$ssc <- ssc(adj = FALSE)

  fit <- function(a) switch(fam,
    poisson = {
      a$fixef.rm <- if (keep_single) "perfect" else "both"
      do.call(fepois, a)
    },
    negbin = {
      a$fixef.rm <- if (keep_single) "perfect" else "both"
      do.call(fenegbin, a)
    },
    {
      a$fixef.rm <- if (keep_single) "none" else "singleton"
      do.call(feols, a)
    })
  if (same_na) {
    M <- fit(args)
    models <- if (inherits(M, "fixest_multi")) as.list(M) else list(M)
  } else {
    models <- lapply(Y, function(v) {
      a <- args
      a$fml <- make_fml(v)
      fit(a)
    })
  }
  t_fit <- proc.time()[["elapsed"]]
  n_d <- nrow(d)
  num <- function(x) format(x, digits = 17L, trim = TRUE, scientific = TRUE)
  w_norm <- function(used) if (length(W)) length(used) / sum(d[[W]][used]) else 1

  for (k in seq_along(models)) {
    m <- models[[k]]
    V <- vcov(m)
    b <- coef(m)
    N <- nobs(m)
    K <- m$nparams
    nb <- length(b)

    pos <- c(X, if (!length(FE)) "(Intercept)")
    bf <- setNames(numeric(length(pos)), pos)
    Vf <- matrix(0, length(pos), length(pos), dimnames = list(pos, pos))
    est <- intersect(pos, names(b))
    bf[est] <- b[est]
    Vf[est, est] <- V[est, est]

    used <- obs(m)
    G <- vapply(CL, function(v) length(unique(d[[v]][used])), 0)

    if (fam == "ols") {
      r <- r2(m, c("r2", "ar2", "wr2", "war2"))
      stat <- c(N = N, K = K, df_a = K - nb, ssr = m$ssr,
                r2 = r[["r2"]], r2_a = r[["ar2"]], r2_within = r[["wr2"]], r2_a_within = r[["war2"]],
                rmse = sqrt(m$ssr * w_norm(used) / (N - K)), ll = NA_real_, r2_p = NA_real_)
    } else {
      r <- r2(m, "pr2")
      stat <- c(N = N, K = K, df_a = K - nb, ssr = NA_real_,
                r2 = NA_real_, r2_a = NA_real_, r2_within = NA_real_, r2_a_within = NA_real_,
                rmse = NA_real_, ll = as.numeric(logLik(m)), r2_p = r[["pr2"]])
    }
    if (length(G)) stat <- c(stat, setNames(G, paste0("G", seq_along(G))))
    stat <- c(stat, n_collin = length(pos) - length(est))

    writeLines(c(
      paste("s", names(stat), num(unname(stat))),
      paste("b", num(unname(bf))),
      paste("v", num(as.vector(Vf)))
    ), file.path(dir, paste0("r", k, ".txt")))

    fesave <- meta$fesave == "1"
    doresid <- one("resid") == "1"
    if (k == 1L && (any(fesave) || doresid)) {
      cols <- list(as.double(d[[one("rid")]][used]))
      if (doresid) cols <- c(cols, list(as.double(resid(m))))
      if (any(fesave)) {
        fv <- fixef(m, notes = FALSE, sorted = FALSE)
        for (j in which(fesave)) cols <- c(cols, list(as.double(fv[[j]][m$fixef_id[[j]]])))
      }
      con <- file(file.path(dir, "save1.bin"), "wb")
      writeBin(unlist(cols), con, size = 8L, endian = "little")
      close(con)
      writeLines(as.character(c(length(used), length(cols))), file.path(dir, "save1.txt"))
    }

    ex <- d[[one("rid")]][setdiff(seq_len(n_d), used)]
    if (length(ex)) writeLines(as.character(ex), file.path(dir, paste0("x", k, ".txt")))
  }

  t_end <- proc.time()[["elapsed"]]
  writeLines(paste(c("lib", "read", "fit", "write"),
                   c(t_lib - t_start, t_read - t_lib, t_fit - t_read, t_end - t_fit)),
             file.path(dir, "time.txt"))
  "ok"
}, error = function(e) conditionMessage(e))

if (!identical(res, "ok")) fail(paste("Error in R:", res))
writeLines("ok", file.path(dir, "ok.txt"))
sink(type = "message")
sink()
close(logcon)
