# Usage:
# Rscript QC_make_html.R <Rscripts directory> <metagStats.txt filepath> <output filepath>

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("The R package 'rmarkdown' is required.", call. = FALSE)
}

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 3L) {
    stop(
        "Usage: QC_make_html.R <R scripts directory> <metagStats.txt> <output HTML>",
        call. = FALSE
    )
}

# Resolve paths before changing context. Previously, relative input/output paths
# were silently reinterpreted relative to the R scripts directory.
start_dir <- getwd()
absolute_path <- function(path) {
    if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) {
        path <- file.path(start_dir, path)
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

scripts_dir <- absolute_path(argv[[1L]])
stats_path <- absolute_path(argv[[2L]])
output_path <- absolute_path(argv[[3L]])
template_path <- file.path(scripts_dir, "QC_html_report.Rmd")

if (!dir.exists(scripts_dir)) {
    stop(sprintf("R scripts directory does not exist: %s", scripts_dir), call. = FALSE)
}
if (!file.exists(template_path)) {
    stop(sprintf("QC report template does not exist: %s", template_path), call. = FALSE)
}
if (!file.exists(stats_path)) {
    stop(sprintf("QC statistics file does not exist: %s", stats_path), call. = FALSE)
}

output_dir <- dirname(output_path)
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop(sprintf("Could not create output directory: %s", output_dir), call. = FALSE)
}

rmarkdown::render(
    input = template_path,
    params = list(statspath = stats_path),
    output_file = basename(output_path),
    output_dir = output_dir,
    knit_root_dir = scripts_dir,
    envir = new.env(parent = globalenv())
)
