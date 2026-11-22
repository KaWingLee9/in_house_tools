# Rebuild the RNA assay, cell metadata, identities, and reductions in Seurat v4.
# seurat_v5: Seurat object loaded with readRDS(), with in-memory RNA matrices.
# Returns a new Seurat object built with the current SeuratObject constructors.
V5toV4 = function(seurat_v5) {
  if (class(seurat_v5)[1L] != "Seurat") {
    stop("seurat_v5 must be a Seurat object.")
  }

  # Read serialized slots directly because Assay5 and LogMap are absent in v4.
  rna = attr(seurat_v5, "assays", exact = TRUE)[["RNA"]]
  if (is.null(rna)) {
    stop("seurat_v5 must contain an RNA assay.")
  }
  meta = attr(seurat_v5, "meta.data", exact = TRUE)
  cells = rownames(meta)
  assay_class = class(rna)[1L]

  if (assay_class %in% c("Assay5", "Assay5T")) {
    layers = attr(rna, "layers", exact = TRUE)
    feature_map = unclass(attr(rna, "features", exact = TRUE))
    cell_map = unclass(attr(rna, "cells", exact = TRUE))

    read_layer = function(layer) {
      layer_names = names(layers)
      parts = layer_names[layer_names == layer | startsWith(layer_names, paste0(layer, "."))]
      if (length(parts) == 0L) {
        return(NULL)
      }
      if (!all(parts %in% colnames(feature_map)) || !all(parts %in% colnames(cell_map))) {
        stop("Missing cell or feature mappings for RNA ", layer, " layers.")
      }
      features = rownames(feature_map)[rowSums(feature_map[, parts, drop = FALSE]) > 0L]
      parts = lapply(parts, function(part) {
        mat = layers[[part]]
        if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
          stop("RNA layer ", part, " must be an in-memory matrix. Materialize it before saving the v5 RDS.")
        }
        if (assay_class == "Assay5T") {
          mat = t(mat)
        }
        layer_features = rownames(feature_map)[feature_map[, part]]
        layer_cells = rownames(cell_map)[cell_map[, part]]
        if (!identical(dim(mat), c(length(layer_features), length(layer_cells)))) {
          stop("Dimensions disagree with the cell and feature mappings for RNA layer ", part, ".")
        }
        dimnames(mat) = list(layer_features, layer_cells)
        if (!identical(layer_features, features)) {
          # Align split matrices by feature ID, inserting sparse zeros where absent.
          alignment = Matrix::sparseMatrix(
            i = match(layer_features, features), j = seq_along(layer_features), x = 1,
            dims = c(length(features), length(layer_features)), dimnames = list(features, layer_features)
          )
          mat = alignment %*% mat
        }
        mat
      })
      layer_cells = unlist(lapply(parts, colnames), use.names = FALSE)
      if (anyDuplicated(layer_cells)) {
        stop("RNA ", layer, " layers contain overlapping cells.")
      }
      do.call(cbind, parts)
    }
    counts = read_layer("counts")
    normalized = read_layer("data")
  } else if (assay_class == "Assay") {
    counts = attr(rna, "counts", exact = TRUE)
    normalized = attr(rna, "data", exact = TRUE)
    if (length(counts) == 0L) counts = NULL
    if (length(normalized) == 0L) normalized = NULL
  } else {
    stop("Unsupported RNA assay class: ", assay_class)
  }

  if (is.null(counts) && is.null(normalized)) {
    stop("The RNA assay must contain a counts or data layer.")
  }
  if ((!is.null(counts) && !setequal(cells, colnames(counts))) ||
      (!is.null(normalized) && !setequal(cells, colnames(normalized)))) {
    stop("Each RNA expression layer must cover the cells in the metadata.")
  }

  if (is.null(counts)) {
    rna_v4 = SeuratObject::CreateAssayObject(data = normalized[, cells, drop = FALSE])
  } else {
    rna_v4 = SeuratObject::CreateAssayObject(counts = counts[, cells, drop = FALSE])
    if (!is.null(normalized)) {
      if (!setequal(rownames(counts), rownames(normalized))) {
        stop("RNA counts and data must contain the same genes to share a v4 assay.")
      }
      rna_v4 = SeuratObject::SetAssayData(
        rna_v4, slot = "data", new.data = normalized[rownames(counts), cells, drop = FALSE]
      )
    }
  }

  seurat_v4 = SeuratObject::CreateSeuratObject(
    counts = rna_v4, assay = "RNA", meta.data = meta,
    project = attr(seurat_v5, "project.name", exact = TRUE)
  )

  identities = attr(seurat_v5, "active.ident", exact = TRUE)
  if (!is.factor(identities) || anyDuplicated(names(identities)) ||
      !setequal(names(identities), cells)) {
    stop("The original Idents must contain one named identity per cell.")
  }
  seurat_v4@active.ident = identities[cells]

  reductions = attr(seurat_v5, "reductions", exact = TRUE)
  if (length(reductions) > 0L) {
    reduction_names = names(reductions)
    if (length(reduction_names) != length(reductions) || anyNA(reduction_names) ||
        any(!nzchar(reduction_names)) || anyDuplicated(reduction_names)) {
      stop("The original reductions must have unique names.")
    }

    # Insert the list directly so reductions computed on a subset keep only those cells.
    seurat_v4@reductions = stats::setNames(lapply(reduction_names, function(reduction_name) {
      reduction = reductions[[reduction_name]]
      embeddings = attr(reduction, "cell.embeddings", exact = TRUE)
      reduction_cells = rownames(embeddings)
      if (!inherits(reduction, "DimReduc") || !is.matrix(embeddings) ||
          length(reduction_cells) == 0L || anyDuplicated(reduction_cells) ||
          !all(reduction_cells %in% cells)) {
        stop("Invalid embeddings or cell IDs in reduction ", reduction_name, ".")
      }
      assay = attr(reduction, "assay.used", exact = TRUE)
      # A global reduction can retain its source assay name even when only RNA is rebuilt.
      global = isTRUE(attr(reduction, "global", exact = TRUE)) || !assay %in% "RNA"
      SeuratObject::CreateDimReducObject(
        embeddings = embeddings[cells[cells %in% reduction_cells], , drop = FALSE],
        loadings = attr(reduction, "feature.loadings", exact = TRUE),
        projected = attr(reduction, "feature.loadings.projected", exact = TRUE),
        assay = assay,
        stdev = attr(reduction, "stdev", exact = TRUE),
        key = attr(reduction, "key", exact = TRUE),
        global = global,
        jackstraw = attr(reduction, "jackstraw", exact = TRUE),
        misc = attr(reduction, "misc", exact = TRUE)
      )
    }), reduction_names)
  }
  methods::validObject(seurat_v4, complete = TRUE)
  seurat_v4
}