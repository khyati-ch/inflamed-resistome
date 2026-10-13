# Exploratory TBK1 checks (not pre-specified)
x2 <- x |> mutate(bin = ntile(log10(baseMean + 1), 10))
b  <- x2$bin[x2$`...1` == 'TBK1']
mean(abs(x2$stat[x2$bin == b]) < abs(x2$stat[x2$`...1` == 'TBK1']))   # 0.926
mean(abs(x$stat) < abs(x$stat[x$`...1` == 'TBK1']))                   # 0.920
ifn_genes <- unique(rx$gene_symbol[rx$gs_name == 'REACTOME_INTERFERON_SIGNALING'])
bg_uc <- uc$stats[!names(uc$stats) %in% ifn_genes]
mean(abs(bg_uc) < abs(uc$stats['TBK1']))                              # 0.921
