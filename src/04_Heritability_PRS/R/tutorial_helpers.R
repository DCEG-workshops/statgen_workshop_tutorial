# Supporting functions for 04_Heritability_PRS.Rmd.

required_inputs <- function(data_dir, phenotype_path) {
  bed_files <- function(prefix) paste0(file.path(data_dir, prefix), c(".bed", ".bim", ".fam"))
  chromosome_files <- function(prefix, suffix) {
    unlist(lapply(1:22, function(chr) paste0(file.path(data_dir, prefix), chr, suffix)))
  }
  list(
    gcta = c(bed_files("chr22"), phenotype_path),
    ldsc = c(file.path(data_dir, c("overall_bc", "lua_bc", "tn_bc")),
             file.path(data_dir, "eur_w_ld_chr/w_hm3.snplist"),
             chromosome_files("eur_w_ld_chr/", c(".l2.ldscore.gz", ".l2.M_5_50"))),
    sldsc = c(chromosome_files("1000G_Phase3_baselineLD_ldscores/baselineLD.",
                             c(".l2.ldscore.gz", ".l2.M_5_50", ".annot.gz")),
              chromosome_files("1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.", ".l2.ldscore.gz"),
              chromosome_files("1000G_Phase3_frq/1000G.EUR.QC.", ".frq")),
    prs = c(bed_files("1kg_eur_22/chr_22"), bed_files("prs_genotype/chr22_test"),
            file.path(data_dir, c("EUR_sum_data", "y_out")))
  )
}

inspect_inputs <- function(requirements) {
  result <- do.call(rbind, lapply(names(requirements), function(section) {
    paths <- requirements[[section]]
    info <- file.info(paths)
    data.frame(section = section, path = paths, bytes = info$size,
               readable = !is.na(info$size) & !info$isdir &
                 info$size > 0 & file.access(paths, 4) == 0,
               stringsAsFactors = FALSE)
  }))
  rownames(result) <- NULL
  result
}

resolve_executable <- function(command) {
  if (!nzchar(command)) return("")
  if (grepl("/", command, fixed = TRUE)) {
    candidate <- path.expand(command)
    if (file.exists(candidate) && file.access(candidate, 1) == 0) {
      return(normalizePath(candidate))
    }
    return("")
  }
  unname(Sys.which(command))
}

command_parts <- function(executable, args) {
  if (endsWith(executable, ".sh")) {
    list(executable = "/bin/bash", args = c(executable, args))
  } else {
    list(executable = executable, args = args)
  }
}

run_command <- function(executable, args, log_prefix, expected = character()) {
  if (!is.character(args) || anyNA(args)) stop("Command arguments must be character values without NA.")
  dir.create(dirname(log_prefix), recursive = TRUE, showWarnings = FALSE)
  parts <- command_parts(executable, args)
  command_text <- paste(shQuote(parts$executable), paste(shQuote(parts$args), collapse = " "))
  writeLines(command_text, paste0(log_prefix, ".command.txt"))
  stdout_file <- paste0(log_prefix, ".stdout.txt")
  stderr_file <- paste0(log_prefix, ".stderr.txt")
  cat("Running", basename(log_prefix), "\n")
  started <- Sys.time()
  status <- system2(parts$executable, shQuote(parts$args), stdout = stdout_file, stderr = stderr_file)
  writeLines(c(paste("started:", started), paste("ended:", Sys.time()),
               paste("exit_status:", status)), paste0(log_prefix, ".status.txt"))
  if (status != 0L) {
    details <- tail(readLines(stderr_file, warn = FALSE), 20L)
    stop("Command failed (exit ", status, "). See ", log_prefix,
         ".stdout.txt and .stderr.txt.\n", paste(details, collapse = "\n"))
  }
  sizes <- file.info(expected)$size
  if (length(expected) && any(is.na(sizes) | sizes <= 0)) {
    stop("Command returned success but an expected output is missing or empty: ",
         paste(expected[is.na(sizes) | sizes <= 0], collapse = ", "))
  }
  invisible(status)
}

read_ldsc_log <- function(prefix) {
  cat(paste(readLines(paste0(prefix, ".log"), warn = FALSE), collapse = "\n"), "\n")
}

validate_gcta_phenotype <- function(fam_path, phenotype_path) {
  fam <- data.table::fread(fam_path, header = FALSE, colClasses = "character")
  pheno <- data.table::fread(phenotype_path, header = FALSE, colClasses = "character")
  if (ncol(fam) != 6L || ncol(pheno) != 3L) {
    stop("GCTA requires a six-column FAM and a headerless three-column FID IID phenotype file.")
  }
  fam_ids <- paste(fam[[1]], fam[[2]], sep = "\t")
  pheno_ids <- paste(pheno[[1]], pheno[[2]], sep = "\t")
  if (anyDuplicated(fam_ids) || anyDuplicated(pheno_ids)) stop("Duplicate GCTA FID/IID pairs.")
  if (!setequal(fam_ids, pheno_ids)) stop("GCTA phenotype FID/IID pairs do not exactly match chr22.fam.")
  values <- suppressWarnings(as.numeric(pheno[[3]]))
  if (any(!is.finite(values)) || any(values == -9) || stats::var(values) <= 0) {
    stop("GCTA phenotypes must be finite, nonmissing, and variable; -9 denotes missing.")
  }
  data.frame(samples = length(values), mean = mean(values), sd = stats::sd(values))
}

initialize_tutorial <- function(params, project_dir) {
  execute <- isTRUE(params$execute)
  if (execute && !nzchar(Sys.getenv("SLURM_JOB_ID"))) {
    stop("Run the analyses inside the Biowulf OnDemand RStudio allocation. ",
         "For a code-only preview, render with params = list(execute = FALSE).")
  }
  data_dir <- normalizePath(path.expand(params$data_dir), mustWork = FALSE)
  phenotype_path <- path.expand(params$gcta_phenotype)
  if (!nzchar(phenotype_path)) {
    phenotype_path <- file.path(dirname(data_dir), "result", "phenotype.phen")
  }
  phenotype_path <- normalizePath(phenotype_path, mustWork = FALSE)
  requested <- c(gcta = params$run_gcta, ldsc = params$run_ldsc,
                 sldsc = params$run_sldsc, prs = params$run_prs)
  if (requested[["sldsc"]] && !requested[["ldsc"]]) {
    stop("run_sldsc = TRUE requires run_ldsc = TRUE to prepare the overall summary statistics.")
  }
  threads <- as.integer(params$threads)
  allocated <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  if (is.na(threads) || threads < 1L || (execute && !is.na(allocated) && threads > allocated)) {
    stop("threads must be positive and no larger than SLURM_CPUS_PER_TASK.")
  }
  data.table::setDTthreads(threads)
  Sys.setenv(OMP_NUM_THREADS = threads, OPENBLAS_NUM_THREADS = threads, MKL_NUM_THREADS = threads)
  output_root <- path.expand(params$output_dir)
  if (!grepl("^/", output_root)) output_root <- file.path(project_dir, output_root)
  output_root <- normalizePath(output_root, mustWork = FALSE)
  shared_root <- normalizePath(dirname(data_dir), mustWork = FALSE)
  if (identical(output_root, shared_root) || startsWith(output_root, paste0(shared_root, "/"))) {
    stop("Choose an output directory outside the shared workshop input directory.")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  run_dir <- tempfile(paste0(if (execute) "run-" else "preview-",
                             format(Sys.time(), "%Y%m%d-%H%M%S"), "-"), tmpdir = output_root)
  if (!dir.create(run_dir)) stop("Could not create a new results directory: ", run_dir)
  requirements <- required_inputs(data_dir, phenotype_path)
  inventory <- inspect_inputs(requirements)
  utils::write.csv(inventory, file.path(run_dir, "input_inventory.csv"), row.names = FALSE)

  if (!params$tool_mode %in% c("biowulf_modules", "path")) stop("Unknown tool_mode.")
  tool_paths <- c(gcta = params$gcta_bin, plink = params$plink_bin, plink2 = params$plink2_bin)
  ldsc_dir <- path.expand(params$ldsc_dir)
  python <- params$ldsc_python
  if (params$tool_mode == "biowulf_modules") {
    tool_paths <- setNames(file.path(project_dir, "scripts", paste0(names(tool_paths), ".sh")), names(tool_paths))
    ldsc_launcher <- file.path(project_dir, "scripts", "ldsc.sh")
  } else {
    tool_paths <- vapply(tool_paths, resolve_executable, character(1))
    python <- resolve_executable(python)
    ldsc_launcher <- python
  }
  tool_ok <- setNames(rep(FALSE, 4L), c("gcta", "plink", "plink2", "ldsc"))
  if (execute) {
    needed <- c(gcta = requested[["gcta"]], plink = requested[["prs"]],
                plink2 = requested[["prs"]], ldsc = requested[["ldsc"]])
    for (name in names(tool_ok)[needed]) {
      if (params$tool_mode == "biowulf_modules") {
        executable <- if (name == "ldsc") ldsc_launcher else tool_paths[[name]]
        parts <- command_parts(executable, "--check")
        status <- system2(parts$executable, shQuote(parts$args),
                          stdout = file.path(run_dir, paste0("tool-", name, ".txt")),
                          stderr = file.path(run_dir, paste0("tool-", name, ".stderr.txt")))
        tool_ok[[name]] <- identical(status, 0L)
      } else {
        tool_ok[[name]] <- if (name == "ldsc") {
          nzchar(python) && nzchar(ldsc_dir) &&
            all(file.exists(file.path(ldsc_dir, c("ldsc.py", "munge_sumstats.py"))))
        } else nzchar(tool_paths[[name]])
      }
    }
  }
  section_tools <- c(gcta = tool_ok[["gcta"]], ldsc = tool_ok[["ldsc"]],
                     sldsc = tool_ok[["ldsc"]], prs = all(tool_ok[c("plink", "plink2")]))
  ready <- vapply(names(requested), function(s) all(inventory$readable[inventory$section == s]), logical(1))
  can_run <- execute & requested & ready & section_tools
  can_run[["sldsc"]] <- can_run[["sldsc"]] && can_run[["ldsc"]]
  status <- data.frame(section = names(requested), requested = unname(requested),
                       inputs_readable = unname(ready), tools_available = unname(section_tools),
                       status = ifelse(can_run, "ready", "not run"))
  utils::write.csv(status, file.path(run_dir, "section_status.csv"), row.names = FALSE)
  saveRDS(params, file.path(run_dir, "parameters.rds"))
  writeLines(c(paste("host:", Sys.info()[["nodename"]]),
               paste("Slurm job:", Sys.getenv("SLURM_JOB_ID")),
               paste("module mode:", params$tool_mode),
               paste(names(tool_paths), tool_paths), paste("LDSC directory:", ldsc_dir),
               paste("LDSC Python/launcher:", ldsc_launcher)), file.path(run_dir, "runtime.txt"))
  list(data_dir = data_dir, phenotype_path = phenotype_path, run_dir = run_dir,
       threads = threads, tool_paths = tool_paths, ldsc_dir = ldsc_dir,
       ldsc_launcher = ldsc_launcher, can_run = can_run, status = status, inventory = inventory)
}
