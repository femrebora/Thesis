#!/usr/bin/env Rscript
# Module 2: STRING PPI Network & Hub Gene Identification
# =================================================================
# Uses STRING REST API (v12.0) directly for lightweight queries.
# Falls back gracefully if STRING API is unreachable.

# Resolve the script directory for both Rscript and source() invocation
script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))),
  error = function(e) tryCatch(dirname(parent.frame(2)$ofile), error = function(e2) getwd())
)
if (is.null(script_dir) || script_dir == "" || script_dir == ".") script_dir <- getwd()
setwd(script_dir)
source("config.R")
source("utils.R")

library(STRINGdb)
library(igraph)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(httr)
library(jsonlite)

dir.create(file.path(RESULTS_DIR, "02_ppi_network"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(RESULTS_DIR, "02_ppi_network")

STRING_API <- "https://version-12-0.string-db.org/api"

# ---- 1. Map symbols to STRING IDs via STRINGdb (lightweight, only aliases) ----
message("[1/5] Connecting to STRING database...")
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  input_directory = ""
)

message("[2/5] Mapping genes to STRING identifiers...")
gene_mapping <- string_db$map(
  data.frame(gene = ALL_GENES),
  "gene",
  removeUnmappedRows = TRUE
)

message("  Mapped ", nrow(gene_mapping), "/", length(ALL_GENES), " genes")
unmapped <- setdiff(ALL_GENES, gene_mapping$gene)
if (length(unmapped) > 0) {
  message("  Unmapped: ", paste(unmapped, collapse = ", "))
}

if (nrow(gene_mapping) == 0) {
  stop("No genes mapped to STRING. Check gene symbols.")
}

# ---- 2. Get PPI network from STRING REST API ----
message("[3/5] Querying STRING REST API for PPI network...")

string_ids <- gene_mapping$STRING_id
symbols <- gene_mapping$gene

# Step A: Get direct interactions among seed genes
identifiers <- paste(string_ids, collapse = "%0d")

get_string_json <- function(endpoint, params, timeout_sec = 60) {
  url <- sprintf("%s/%s?%s&caller_identity=gbm_analysis", STRING_API, endpoint, params)
  resp <- GET(url, timeout(timeout_sec))
  if (status_code(resp) != 200) {
    stop(sprintf("STRING API returned %d for %s", status_code(resp), endpoint))
  }
  fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
}

# Step A: Network among seed genes
params_network <- sprintf("identifiers=%s&species=9606&required_score=400", identifiers)
ppi_data <- get_string_json("json/network", params_network)
message("  Direct edges among seeds: ", nrow(ppi_data))

# Step B: Get interaction partners (1st shell) for each seed gene
# Limit to top 10 partners per query to keep it manageable
interactor_data <- data.frame(
  stringId_A = character(), stringId_B = character(),
  preferredName_A = character(), preferredName_B = character(),
  score = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(string_ids)) {
  sid <- string_ids[i]
  sym <- symbols[i]
  params_partners <- sprintf(
    "identifiers=%s&species=9606&required_score=400&limit=10",
    URLencode(sid, reserved = TRUE)
  )
  tryCatch({
    partners <- get_string_json("json/interaction_partners", params_partners)
    if (nrow(partners) > 0) {
      interactor_data <- rbind(interactor_data, partners)
    }
    message("  ", sym, ": ", nrow(partners), " partners")
  }, error = function(e) {
    message("  ", sym, ": API error - ", e$message)
  })
  Sys.sleep(0.3)  # rate limiting courtesy
}

message("  Total interaction partners: ", nrow(interactor_data))

# Combine: use direct edges (preferred) + interaction partners
combined <- rbind(
  data.frame(
    from = ppi_data$preferredName_A,
    to = ppi_data$preferredName_B,
    score = ppi_data$score / 1000,
    stringsAsFactors = FALSE
  ),
  data.frame(
    from = interactor_data$preferredName_A,
    to = interactor_data$preferredName_B,
    score = interactor_data$score / 1000,
    stringsAsFactors = FALSE
  )
)

# Remove self-loops
combined <- combined[combined$from != combined$to, ]

if (nrow(combined) == 0) {
  stop("No PPI edges retrieved from STRING API. The API may be unreachable.")
}

# Build graph, collapsing multi-edges (take max score)
g <- graph_from_data_frame(combined, directed = FALSE)
g <- igraph::simplify(g, edge.attr.comb = list(score = "max", weight = "max"))
E(g)$weight <- E(g)$score

# Identify seed genes
V(g)$symbol <- V(g)$name
V(g)$is_seed <- V(g)$symbol %in% symbols

message("  Final network: ", vcount(g), " nodes, ", ecount(g), " edges")
message("  Seed genes in network: ", sum(V(g)$is_seed))

# ---- 3. Hub gene identification ----
message("[4/5] Computing hub metrics...")

V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)
V(g)$closeness <- closeness(g, normalized = TRUE)
V(g)$eigenvector <- eigen_centrality(g)$vector

# Simplified MCC: degree * (local_clustering + 1)
V(g)$local_clustering <- transitivity(g, type = "local", isolates = "zero")
V(g)$mcc_score <- V(g)$degree * (V(g)$local_clustering + 1)

hub_df <- data.frame(
  symbol = V(g)$symbol,
  is_seed = V(g)$is_seed,
  degree = V(g)$degree,
  betweenness = round(V(g)$betweenness, 4),
  closeness = round(V(g)$closeness, 4),
  eigenvector = round(V(g)$eigenvector, 4),
  mcc_score = round(V(g)$mcc_score, 2),
  stringsAsFactors = FALSE
)
hub_df <- hub_df[order(-hub_df$mcc_score), ]
hub_df$rank <- 1:nrow(hub_df)

message("  Top 5 hub genes: ", paste(head(hub_df$symbol[hub_df$is_seed], 5), collapse = ", "))

# ---- 4. Community detection and GO enrichment ----
message("[5/5] Detecting modules and running enrichment...")

communities <- cluster_walktrap(g, steps = 4)
V(g)$module <- communities$membership

module_sizes <- table(communities$membership)
enrichment_results <- list()

for (mod_id in as.numeric(names(module_sizes[module_sizes >= 5]))) {
  mod_genes <- V(g)$symbol[V(g)$module == mod_id]
  mod_genes <- mod_genes[mod_genes != ""]

  tryCatch({
    ego <- enrichGO(
      gene = mod_genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      ego_df <- as.data.frame(ego)
      ego_df$module <- mod_id
      ego_df$module_size <- length(mod_genes)
      enrichment_results[[length(enrichment_results) + 1]] <- head(ego_df, 5)
    }
  }, error = function(e) {
    message("  Enrichment failed for module ", mod_id, ": ", e$message)
  })
}

if (length(enrichment_results) > 0) {
  module_enrich <- do.call(rbind, enrichment_results)
} else {
  module_enrich <- data.frame(
    module = integer(), Description = character(), p.adjust = numeric(),
    message = "No enriched GO terms found at q < 0.05"
  )
}

# ---- 5. PPI network visualization ----
message("[5/5] Plotting PPI network...")

V(g)$color <- ifelse(V(g)$is_seed, "#E41A1C", "#999999")
V(g)$size <- 2 + log2(V(g)$degree + 1) * 2
top_nodes <- names(sort(V(g)$degree, decreasing = TRUE))[1:min(20, vcount(g))]
V(g)$label <- ifelse(V(g)$name %in% top_nodes, V(g)$symbol, "")
V(g)$label.cex <- 0.6

set.seed(42)
layout <- layout_with_fr(g, niter = 2000)

pdf(file.path(OUT, "ppi_network.pdf"), width = 12, height = 10)
par(mar = c(1, 1, 2, 1))
plot(g,
  layout = layout,
  vertex.color = V(g)$color,
  vertex.size = V(g)$size,
  vertex.label = V(g)$label,
  vertex.label.cex = V(g)$label.cex,
  vertex.frame.color = NA,
  edge.width = E(g)$weight * 3,
  edge.color = rgb(0.7, 0.7, 0.7, 0.4),
  main = "PPI Network: CRISPR Screen Hits + 1st-Shell Interactors"
)
legend("bottomleft",
  legend = c("Seed (CRISPR hit)", "Interactor"),
  pt.bg = c("#E41A1C", "#999999"), pch = 21, pt.cex = 1.5,
  bty = "n"
)
dev.off()

tiff(file.path(OUT, "ppi_network.tiff"), width = 12, height = 10, units = "in", res = 300, compression = "lzw")
par(mar = c(1, 1, 2, 1))
plot(g,
  layout = layout,
  vertex.color = V(g)$color,
  vertex.size = V(g)$size,
  vertex.label = V(g)$label,
  vertex.label.cex = V(g)$label.cex,
  vertex.frame.color = NA,
  edge.width = E(g)$weight * 3,
  edge.color = rgb(0.7, 0.7, 0.7, 0.4),
  main = "PPI Network: CRISPR Screen Hits + 1st-Shell Interactors"
)
legend("bottomleft",
  legend = c("Seed (CRISPR hit)", "Interactor"),
  pt.bg = c("#E41A1C", "#999999"), pch = 21, pt.cex = 1.5,
  bty = "n"
)
dev.off()

# ---- 6. Save tables ----
save_table(hub_df, file.path(OUT, "hub_genes"))
save_table(module_enrich, file.path(OUT, "module_enrichment"))

message("[done] Module 2 complete. Outputs in ", OUT)
