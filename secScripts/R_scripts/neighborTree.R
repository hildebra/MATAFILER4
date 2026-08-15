if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The R package 'ape' is required.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 2L) {
  stop("Usage: neighborTree.R <tree file> <target tip>|--all", call. = FALSE)
}

inF <- args[[1L]]
target <- args[[2L]]

if (!file.exists(inF)) {
  stop(sprintf("Tree file does not exist: %s", inF), call. = FALSE)
}

tree <- tryCatch(
  ape::read.tree(inF),
  error = function(e) stop(sprintf("Could not read tree '%s': %s", inF, conditionMessage(e)), call. = FALSE)
)
if (inherits(tree, "multiPhylo")) {
  stop("The input must contain exactly one tree.", call. = FALSE)
}
if (is.null(tree) || !inherits(tree, "phylo")) {
  stop(sprintf("No valid tree was found in: %s", inF), call. = FALSE)
}
if (anyDuplicated(tree$tip.label)) {
  stop("Tree tip labels must be unique.", call. = FALSE)
}

ordered_neighbors <- function(distances) {
  neighbors <- sort(distances[is.finite(distances) & distances >= 0.01])
  names(neighbors)
}

if (identical(target, "--all")) {
  # Compute the tip-to-tip distances once.  Phase II consumes one row per MGS,
  # avoiding a separate R startup, tree parse, and distance calculation per tip.
  distances <- ape::cophenetic.phylo(tree)
  for (tip in tree$tip.label) {
    neighbors <- ordered_neighbors(distances[tip, ])
    cat(tip, "\t", paste(neighbors, collapse = " "), "\n", sep = "")
  }
  quit(save = "no", status = 0)
}

if (!target %in% tree$tip.label) {
  stop(sprintf("Target tip '%s' is not present in the tree.", target), call. = FALSE)
}

neighbors <- ordered_neighbors(ape::cophenetic.phylo(tree)[target, ])
# Do not index to an arbitrary minimum length: doing so pads short results with
# NA values, which downstream Perl callers interpret as candidate tip names.
cat(paste(neighbors, collapse = " "), "\n", sep = "")
