# XJC
suppressPackageStartupMessages({
  library(data.table)
  library(STRINGdb)
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

out_dir <- "F:/222/练习/data/processed/GSE160611"
res <- fread(file.path(out_dir, "DE_voom_limma_OS_vs_PS_4hr.tsv"))

# -------- parameters --------
fdr_cutoff <- 0.05      # 显著性阈值；边不够可改 0.1
max_genes <- 200        # 节点数上限；图太大可改 100
score_cutoff <- 700     # STRING combined_score: 400/700/900
layout_seed <- 1

# 1) choose genes
res_sig <- res[!is.na(adj.P.Val) & adj.P.Val < fdr_cutoff]
setorder(res_sig, adj.P.Val)
genes <- unique(res_sig$gene)[1:min(max_genes, uniqueN(res_sig$gene))]

if (length(genes) < 10) {
  stop("Too few genes at current fdr_cutoff. Try fdr_cutoff=0.1 or increase max_genes.")
}

gene_df <- data.frame(gene = genes, stringsAsFactors = FALSE)

# 2) init STRINGdb (human 9606)
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 0,
  input_directory = file.path(out_dir, "STRINGdb_cache")
)

# 3) map symbols -> STRING ids
mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = FALSE)

cat("Mapping summary:\n")
cat("  input genes:   ", nrow(gene_df), "\n", sep = "")
cat("  mapped genes:  ", sum(!is.na(mapped$STRING_id)), "\n", sep = "")
cat("  unmapped genes:", sum(is.na(mapped$STRING_id)), "\n", sep = "")

mapped_ok <- mapped[!is.na(mapped$STRING_id), ]
if (nrow(mapped_ok) < 10) stop("Too few mapped genes. Something is wrong with mapping/species/version.")

# 4) get interactions among mapped proteins
ppi <- string_db$get_interactions(mapped_ok$STRING_id)
keep_ids <- unique(mapped_ok$STRING_id)
ppi <- ppi[ppi$from %in% keep_ids & ppi$to %in% keep_ids, ]
ppi <- ppi[ppi$combined_score >= score_cutoff, ]
if (nrow(ppi) == 0) stop("No edges pass score_cutoff. Try score_cutoff=400 or use more genes.")

# 5) build nodes table
nodes <- as.data.table(mapped_ok)[, .(gene, STRING_id)]
nodes <- merge(nodes, res[, .(gene, logFC, adj.P.Val)], by = "gene", all.x = TRUE)

edges <- as.data.table(ppi)
setnames(edges, c("from", "to", "combined_score"), c("STRING_from", "STRING_to", "score"))
edges[, weight := score / 1000]

id2gene <- nodes[, .(STRING_id, gene)]
setkey(id2gene, STRING_id)
edges[, gene_from := id2gene[.(STRING_from)]$gene]
edges[, gene_to   := id2gene[.(STRING_to)]$gene]

fwrite(nodes, file.path(out_dir, "ppi_nodes_STRING_OS_vs_PS_4hr.tsv.gz"), sep = "\t")
fwrite(edges, file.path(out_dir, "ppi_edges_STRING_OS_vs_PS_4hr.tsv.gz"), sep = "\t")

g <- graph_from_data_frame(
  d = edges[, .(from = STRING_from, to = STRING_to, weight, score)],
  directed = FALSE,
  vertices = nodes[, .(name = STRING_id, gene, logFC, adj.P.Val)]
)

set.seed(layout_seed)
p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.25, color = "grey50") +
  geom_node_point(aes(color = logFC, size = -log10(pmax(adj.P.Val, 1e-300))), alpha = 0.9) +
  scale_color_gradient2(low = "#2c7bb6", mid = "grey90", high = "#d7191c", midpoint = 0) +
  geom_node_text(aes(label = gene), repel = TRUE, size = 3, max.overlaps = 80) +
  theme_void()

ggsave(file.path(out_dir, "ppi_network_STRING_OS_vs_PS_4hr.pdf"), p, width = 10, height = 8)
