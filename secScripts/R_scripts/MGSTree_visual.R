# Plot an MGS phylogeny and automatically annotate phylum-level clades.

suppressPackageStartupMessages({
  library("ape")
  library("ggplot2")
  library("ggtree")
  library("phytools")
})

stop_with_usage <- function() {
  stop(
    "Usage: MGSTree_visual.R <MGS.matL7.txt> <tree file> <output PDF>",
    call. = FALSE
  )
}

read_taxonomy <- function(path, tip_labels) {
  raw <- tryCatch(
    read.delim(
      path,
      header = TRUE,
      row.names = 1L,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) {
      stop(sprintf("Could not read metadata '%s': %s", path, conditionMessage(e)), call. = FALSE)
    }
  )

  taxonomy_names <- c(
    "superkingdom", "phylum", "class", "order",
    "family", "genus", "species", "MGS"
  )
  fields <- strsplit(rownames(raw), ";", fixed = TRUE)
  valid <- lengths(fields) == length(taxonomy_names)
  if (!all(valid)) {
    warning(sprintf(
      "Ignoring %d metadata row(s) that do not contain exactly eight semicolon-separated taxonomy fields.",
      sum(!valid)
    ), call. = FALSE)
    fields <- fields[valid]
  }
  if (!length(fields)) {
    stop("The metadata contains no valid eight-rank taxonomy rows.", call. = FALSE)
  }

  metadata <- as.data.frame(
    do.call(rbind, fields),
    stringsAsFactors = FALSE
  )
  names(metadata) <- taxonomy_names
  metadata[] <- lapply(metadata, trimws)

  metadata <- metadata[grepl("^MGS", metadata$MGS) & nzchar(metadata$phylum), , drop = FALSE]
  if (!nrow(metadata)) {
    stop("The metadata contains no MGS rows with a phylum assignment.", call. = FALSE)
  }

  duplicated_mgs <- duplicated(metadata$MGS)
  if (any(duplicated_mgs)) {
    warning(sprintf(
      "Ignoring %d duplicate MGS metadata row(s); the first assignment is used.",
      sum(duplicated_mgs)
    ), call. = FALSE)
    metadata <- metadata[!duplicated_mgs, , drop = FALSE]
  }

  metadata <- metadata[metadata$MGS %in% tip_labels, , drop = FALSE]
  if (!nrow(metadata)) {
    stop("None of the MGS identifiers in the metadata match a tree tip.", call. = FALSE)
  }

  missing_metadata <- setdiff(tip_labels, metadata$MGS)
  if (length(missing_metadata)) {
    warning(sprintf(
      "%d of %d tree tips have no matching MGS metadata and will remain unannotated.",
      length(missing_metadata), length(tip_labels)
    ), call. = FALSE)
  }

  metadata$species[metadata$species %in% c("", "?")] <- paste(
    metadata$genus[metadata$species %in% c("", "?")],
    "unclass"
  )
  # ggtree's `%<+%` joins metadata onto the tree by renaming the *first*
  # column of this data frame to "label" and matching it against the tree's
  # tip labels. The identifier column (MGS) must therefore come first; and
  # since that rename creates a "label" column itself, we must not also keep
  # a separately-named "label" column, or dplyr::rename() fails with
  # "Names must be unique" (duplicate "label" columns).
  metadata <- metadata[, c("MGS", setdiff(names(metadata), "MGS"))]
  metadata
}

root_tree <- function(tree, metadata) {
  # Group tip labels by superkingdom (e.g. Archaea vs Bacteria) so we can try
  # rooting the tree between whichever superkingdoms are actually present,
  # rather than assuming a single fixed outgroup taxon.
  by_superkingdom <- split(metadata$MGS, metadata$superkingdom)
  by_superkingdom <- lapply(by_superkingdom, intersect, tree$tip.label)
  by_superkingdom <- by_superkingdom[lengths(by_superkingdom) > 0L]

  if (length(by_superkingdom) >= 2L) {
    message(sprintf(
      "Detected %d superkingdom(s) among the tips: %s.",
      length(by_superkingdom), paste(names(by_superkingdom), collapse = ", ")
    ))

    # Try the smallest superkingdom group as the outgroup first (most likely
    # to be genuinely monophyletic and gives the best-supported root
    # placement between the represented superkingdoms); if it isn't
    # monophyletic or rooting otherwise fails, fall through to the next
    # smallest superkingdom.
    ordered <- by_superkingdom[order(lengths(by_superkingdom))]
    for (name in names(ordered)) {
      outgroup <- ordered[[name]]
      if (length(outgroup) >= length(tree$tip.label)) {
        next
      }

      is_mono <- length(outgroup) == 1L || isTRUE(ape::is.monophyletic(tree, outgroup))
      if (!is_mono) {
        warning(sprintf(
          "%s tips are not monophyletic; trying another superkingdom split.",
          name
        ), call. = FALSE)
        next
      }

      rooted <- tryCatch(
        ape::root(tree, outgroup = outgroup, resolve.root = TRUE),
        error = function(e) {
          warning(sprintf(
            "Could not root between superkingdoms using %s as the outgroup (%s).",
            name, conditionMessage(e)
          ), call. = FALSE)
          NULL
        }
      )
      if (!is.null(rooted)) {
        message(sprintf(
          "Tree rooted between superkingdoms using %d %s tip(s) as the outgroup.",
          length(outgroup), name
        ))
        return(rooted)
      }
    }

    warning(
      "Could not root between the represented superkingdoms; using midpoint rooting instead.",
      call. = FALSE
    )
  }

  message("Tree rooted using midpoint rooting.")
  phytools::midpoint.root(tree)
}

# Return maximal subtrees whose annotated descendant tips all have the same
# phylum. This naturally splits paraphyletic/polyphyletic phyla and avoids the
# overlapping-MRCA loops that previously failed on ordinary tree topologies.
find_pure_clades <- function(tree, metadata) {
  tip_count <- length(tree$tip.label)
  node_count <- tip_count + tree$Nnode
  purity <- rep(NA_character_, node_count)
  purity[match(metadata$MGS, tree$tip.label)] <- metadata$phylum

  children <- split(tree$edge[, 2L], tree$edge[, 1L])
  parent <- rep(NA_integer_, node_count)
  parent[tree$edge[, 2L]] <- tree$edge[, 1L]
  remaining <- integer(node_count)
  parent_nodes <- as.integer(names(children))
  remaining[parent_nodes] <- lengths(children)

  # Process tips upward. An internal node is queued only after all of its
  # children have been evaluated, so this works for binary and multifurcating trees.
  queue <- integer(node_count)
  queue[seq_len(tip_count)] <- seq_len(tip_count)
  queue_length <- tip_count
  cursor <- 1L
  while (cursor <= queue_length) {
    node <- queue[[cursor]]
    cursor <- cursor + 1L
    ancestor <- parent[[node]]
    if (is.na(ancestor)) {
      next
    }

    remaining[[ancestor]] <- remaining[[ancestor]] - 1L
    if (remaining[[ancestor]] == 0L) {
      child_values <- purity[children[[as.character(ancestor)]]]
      if (all(!is.na(child_values)) && length(unique(child_values)) == 1L) {
        purity[[ancestor]] <- child_values[[1L]]
      }
      queue_length <- queue_length + 1L
      queue[[queue_length]] <- ancestor
    }
  }

  pure_nodes <- which(!is.na(purity))
  maximal <- vapply(
    pure_nodes,
    function(node) {
      ancestor <- parent[[node]]
      is.na(ancestor) || is.na(purity[[ancestor]]) || purity[[ancestor]] != purity[[node]]
    },
    logical(1L)
  )

  data.frame(
    clade = unname(purity[pure_nodes[maximal]]),
    node = pure_nodes[maximal],
    stringsAsFactors = FALSE
  )
}

main <- function(args) {
  if (length(args) < 3L) {
    stop_with_usage()
  }

  in_meta <- args[[1L]]
  in_tree <- args[[2L]]
  out_pdf <- args[[3L]]
  if (!file.exists(in_meta)) {
    stop(sprintf("Metadata file does not exist: %s", in_meta), call. = FALSE)
  }
  if (!file.exists(in_tree)) {
    stop(sprintf("Tree file does not exist: %s", in_tree), call. = FALSE)
  }

  tree <- tryCatch(
    ape::read.tree(in_tree),
    error = function(e) stop(sprintf("Could not read tree '%s': %s", in_tree, conditionMessage(e)), call. = FALSE)
  )
  if (inherits(tree, "multiPhylo")) {
    stop("The input must contain exactly one tree.", call. = FALSE)
  }
  if (is.null(tree) || !inherits(tree, "phylo") || length(tree$tip.label) < 2L) {
    stop("The input does not contain a valid tree with at least two tips.", call. = FALSE)
  }
  if (anyDuplicated(tree$tip.label)) {
    stop("Tree tip labels must be unique.", call. = FALSE)
  }

  metadata <- read_taxonomy(in_meta, tree$tip.label)
  rooted_tree <- root_tree(tree, metadata)
  clades <- find_pure_clades(rooted_tree, metadata)
  if (!nrow(clades)) {
    stop("No annotated clades could be identified in the tree.", call. = FALSE)
  }

  tree_plot <- ggtree(rooted_tree, layout = "circular")
  tree_xmax <- suppressWarnings(max(tree_plot$data$x, na.rm = TRUE))
  if (!is.finite(tree_xmax) || tree_xmax <= 0) {
    tree_xmax <- 1
  }

  plot0 <- tree_plot %<+% metadata +
    geom_tree(aes(color = phylum), linewidth = 0.8) +
    geom_tiplab(aes(label = label), size = 1.4) +
    xlim(NA, tree_xmax * 1.55) +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(),
      legend.key = element_blank(),
      legend.key.size = grid::unit(0.4, "cm"),
      legend.text = element_text(size = 6),
      title = element_text(size = 8)
    )

  plot1 <- tree_plot %<+% metadata +
    geom_highlight(
      data = clades,
      aes(node = node, fill = clade),
      show.legend = TRUE
    ) +
    geom_tiplab(aes(label = species), size = 1.8) +
    xlim(NA, tree_xmax * 1.55) +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(),
      legend.key = element_blank(),
      legend.key.size = grid::unit(0.4, "cm"),
      legend.text = element_text(size = 5),
      title = element_text(size = 8)
    )

  clade_names <- sort(unique(clades$clade))
  grey_values <- grDevices::gray.colors(
    length(clade_names), start = 0.97, end = 0.72
  )
  names(grey_values) <- clade_names
  label_offset <- tree_xmax * 0.08

  plot2 <- tree_plot %<+% metadata +
    geom_highlight(
      data = clades,
      aes(node = node, fill = clade),
      alpha = 1,
      align = "both",
      extend = tree_xmax * 0.02,
      show.legend = FALSE
    ) +
    geom_cladelab(
      data = clades,
      mapping = aes(node = node, label = clade),
      fontsize = 3,
      align = TRUE,
      angle = "auto",
      offset.text = label_offset
    ) +
    geom_tree(linewidth = 0.3) +
    geom_tippoint() +
    xlim(NA, tree_xmax * 2.1) +
    scale_fill_manual(values = grey_values)

  output_dir <- dirname(out_pdf)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop(sprintf("Could not create output directory: %s", output_dir), call. = FALSE)
  }

  grDevices::pdf(out_pdf)
  device_number <- grDevices::dev.cur()
  on.exit({
    if (device_number %in% grDevices::dev.list()) {
      grDevices::dev.off(device_number)
    }
  }, add = TRUE)
  print(plot0)
  print(plot1)
  print(plot2)
  grDevices::dev.off(device_number)
  message(sprintf("Wrote phylogeny plots to %s", out_pdf))
}

main(commandArgs(trailingOnly = TRUE))