# Discovery analysis: nsp13 interactors vs host expression in IBD
# Inputs: data/zenodo/Table2.xlsx, data/DE_UC_vs_nonIBD.csv, data/DE_CD_vs_nonIBD.csv
library(readxl); library(tidyverse); library(fgsea); library(msigdbr)
set.seed(2026)

# 1. Interaction sets
ppi <- read_excel('data/zenodo/Table2.xlsx', skip = 1) |>
  transmute(bait = sub('^SARS-CoV2 ', '', Bait), gene = PreyGene) |>
  filter(!is.na(gene)) |> distinct()
bait_sets <- split(ppi$gene, ppi$bait)

# 2. Host results
de_uc <- read_csv('data/DE_UC_vs_nonIBD.csv', show_col_types = FALSE)
de_cd <- read_csv('data/DE_CD_vs_nonIBD.csv', show_col_types = FALSE)

# 3. Reactome interferon sets (fixed in advance)
rx  <- msigdbr(species = 'Homo sapiens', collection = 'C2', subcollection = 'CP:REACTOME')
ifn <- c('REACTOME_INTERFERON_SIGNALING', 'REACTOME_ISG15_ANTIVIRAL_MECHANISM')
ifn_sets <- lapply(setNames(ifn, ifn), function(s) unique(rx$gene_symbol[rx$gs_name == s]))

# 4. Analysis function
run_contrast <- function(raw, label) {
  de <- raw |> transmute(gene = `...1`, stat = stat, bm = baseMean) |>
    filter(!is.na(gene), !is.na(stat)) |> distinct(gene, .keep_all = TRUE)
  stats <- setNames(de$stat, de$gene)
  sets <- lapply(bait_sets, intersect, names(stats))
  fg <- fgsea(sets, stats, minSize = 5, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |> mutate(rank_absNES = rank(-abs(NES))) |> arrange(rank_absNES)
  write_csv(select(fg, -leadingEdge), sprintf('results/proteomics/fgsea_baits_%s.csv', label))
  n13 <- sets$nsp13
  de$bin <- ntile(log10(de$bm + 1), 10)
  k <- as.integer(table(factor(de$bin[match(n13, de$gene)], levels = 1:10)))
  rand <- replicate(1000, unlist(lapply(1:10, function(b)
    sample(setdiff(de$gene[de$bin == b], n13), k[b]))), simplify = FALSE)
  names(rand) <- paste0('r', seq_along(rand))
  fr  <- fgsea(rand, stats, minSize = 5, maxSize = 500, nPermSimple = 1000)
  obs <- fg$NES[fg$pathway == 'nsp13']
  p_rand <- (1 + sum(abs(fr$NES) >= abs(obs), na.rm = TRUE)) / (1 + sum(!is.na(fr$NES)))
  ifn_fg <- fgsea(ifn_sets, stats, minSize = 10, maxSize = 1000, nPermSimple = 10000) |> as_tibble()
  ifn_fg$padj_holm <- p.adjust(ifn_fg$pval, 'holm')
  write_csv(select(ifn_fg, -leadingEdge), sprintf('results/proteomics/ifn_sets_%s.csv', label))
  list(fg = fg, stats = stats, sets = sets, p_rand = p_rand, ifn = ifn_fg)
}

# 5. Run both contrasts
uc <- run_contrast(de_uc, 'UC')   # primary
cd <- run_contrast(de_cd, 'CD')   # secondary

# 6. Figure
ggsave('results/proteomics/nsp13_UC_enrichment.png',
       plotEnrichment(uc$sets$nsp13, uc$stats) + ggtitle('nsp13 interactors in UC vs nonIBD'),
       width = 7, height = 4.5, dpi = 300)
