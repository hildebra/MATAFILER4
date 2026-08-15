if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The R package 'ape' is required.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(paste(
    "Usage: neighborTree.R <tree file> <target tip>|--all",
    "[--preferred FILE | --preferred-tip TIP]",
    "[--preferred-quantile FLOAT] [--preferred-nearest-factor FLOAT]"
  ), call. = FALSE)
}

inF <- args[[1L]]
target <- args[[2L]]
preferredFile <- ""
preferredTip <- ""
preferredQuantile <- 0.50
preferredNearestFactor <- 10.0

optionArgs <- args[-seq_len(2L)]
i <- 1L
while (i <= length(optionArgs)) {
  option <- optionArgs[[i]]
  if (!option %in% c(
    "--preferred", "--preferred-tip", "--preferred-quantile",
    "--preferred-nearest-factor"
  )) {
    stop(sprintf("Unknown option: %s", option), call. = FALSE)
  }
  if (i == length(optionArgs)) {
    stop(sprintf("Missing value for option: %s", option), call. = FALSE)
  }
  value <- optionArgs[[i + 1L]]
  if (identical(option, "--preferred")) {
    preferredFile <- value
  } else if (identical(option, "--preferred-tip")) {
    preferredTip <- value
  } else if (identical(option, "--preferred-quantile")) {
    preferredQuantile <- suppressWarnings(as.numeric(value))
  } else if (identical(option, "--preferred-nearest-factor")) {
    preferredNearestFactor <- suppressWarnings(as.numeric(value))
  }
  i <- i + 2L
}

if (nzchar(preferredFile) && nzchar(preferredTip)) {
  stop("Use only one of --preferred and --preferred-tip.", call. = FALSE)
}
if (!is.finite(preferredQuantile) || preferredQuantile <= 0 || preferredQuantile > 1) {
  stop("--preferred-quantile must be greater than 0 and at most 1.", call. = FALSE)
}
if (!is.finite(preferredNearestFactor) || preferredNearestFactor < 1) {
  stop("--preferred-nearest-factor must be at least 1.", call. = FALSE)
}
if (!file.exists(inF)) {
  stop(sprintf("Tree file does not exist: %s", inF), call. = FALSE)
}
if (nzchar(preferredFile) && !file.exists(preferredFile)) {
  stop(sprintf("Preferred-outgroup file does not exist: %s", preferredFile), call. = FALSE)
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

preferredByTarget <- character()
if (nzchar(preferredFile)) {
  preferredLines <- readLines(preferredFile, warn = FALSE)
  for (lineNumber in seq_along(preferredLines)) {
    line <- sub("[\\r\\n]+$", "", preferredLines[[lineNumber]])
    if (!nzchar(line) || grepl("^\\s*#", line)) {
      next
    }
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) < 2L || !nzchar(fields[[1L]]) || !nzchar(fields[[2L]])) {
      stop(sprintf("Malformed preferred-outgroup row %d in %s", lineNumber, preferredFile), call. = FALSE)
    }
    if (fields[[1L]] %in% names(preferredByTarget)) {
      stop(sprintf("Duplicate preferred-outgroup target '%s' in %s", fields[[1L]], preferredFile), call. = FALSE)
    }
    preferredByTarget[[fields[[1L]]]] <- fields[[2L]]
  }
}

ranked_neighbors <- function(distances, targetTip, preferred = "") {
  eligible <- sort(distances[is.finite(distances) & distances >= 0.01])
  candidateNames <- names(eligible)
  result <- list(
    candidates = candidateNames,
    decision = "none",
    preferred = preferred,
    preferredDistance = NA_real_,
    cutoff = NA_real_
  )
  if (!nzchar(preferred)) {
    return(result)
  }
  if (identical(preferred, targetTip)) {
    result$decision <- "same_as_target"
    return(result)
  }
  if (!preferred %in% names(distances)) {
    result$decision <- "absent_from_tree"
    return(result)
  }
  preferredDistance <- unname(distances[[preferred]])
  result$preferredDistance <- preferredDistance
  if (!is.finite(preferredDistance)) {
    result$decision <- "non_finite_distance"
    return(result)
  }
  if (preferredDistance < 0.01) {
    result$decision <- "too_close"
    return(result)
  }
  if (!length(eligible)) {
    result$decision <- "no_eligible_neighbors"
    return(result)
  }
  nearestDistance <- min(unname(eligible))
  quantileCutoff <- unname(stats::quantile(
    unname(eligible), probs = preferredQuantile, names = FALSE, type = 8
  ))
  cutoff <- max(quantileCutoff, nearestDistance * preferredNearestFactor)
  result$cutoff <- cutoff
  if (preferredDistance > cutoff) {
    result$decision <- "too_distant"
    result$candidates <- candidateNames[candidateNames != preferred]
    return(result)
  }
  result$decision <- "accepted"
  result$candidates <- c(preferred, candidateNames[candidateNames != preferred])
  result
}

format_number <- function(value) {
  if (!is.finite(value)) "" else format(value, digits = 10L, trim = TRUE, scientific = FALSE)
}

if (identical(target, "--all")) {
  # Compute tip-to-tip distances once. Phase II consumes one authoritative row
  # per MGS, avoiding a separate R startup, tree parse, and distance calculation.
  distances <- ape::cophenetic.phylo(tree)
  for (tip in tree$tip.label) {
    preferred <- if (tip %in% names(preferredByTarget)) {
      unname(preferredByTarget[tip])
    } else {
      ""
    }
    ranked <- ranked_neighbors(distances[tip, ], tip, preferred)
    cat(
      tip, "\t", ranked$decision, "\t", ranked$preferred, "\t",
      format_number(ranked$preferredDistance), "\t", format_number(ranked$cutoff), "\t",
      paste(ranked$candidates, collapse = " "), "\n", sep = ""
    )
  }
  quit(save = "no", status = 0)
}

if (!target %in% tree$tip.label) {
  stop(sprintf("Target tip '%s' is not present in the tree.", target), call. = FALSE)
}
if (!nzchar(preferredTip)) {
  preferredTip <- if (target %in% names(preferredByTarget)) {
    unname(preferredByTarget[target])
  } else {
    ""
  }
}
ranked <- ranked_neighbors(ape::cophenetic.phylo(tree)[target, ], target, preferredTip)
# Single-target mode remains backward compatible: callers receive only the
# authoritative, ordered candidate list.
cat(paste(ranked$candidates, collapse = " "), "\n", sep = "")
