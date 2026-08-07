if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The R package 'ape' is required.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: neighborTree.R <tree file> <target tip>", call. = FALSE)
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
if (!target %in% tree$tip.label) {
  stop(sprintf("Target tip '%s' is not present in the tree.", target), call. = FALSE)
}

distances <- ape::cophenetic.phylo(tree)[target, ]
neighbors <- sort(distances[is.finite(distances) & distances >= 0.01])

# Do not index to an arbitrary minimum length: doing so pads short results with
# NA values, which downstream Perl callers interpret as candidate tip names.
cat(paste(names(neighbors), collapse = " "), "\n", sep = "")





