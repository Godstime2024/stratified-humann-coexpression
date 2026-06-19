# =============================================================================
# h2_sink_compare.R
# Where does H2 go when methanogenesis is suppressed?
#
# Tests CON vs SCFP for alternative SCFA / electron sink pathways:
#   - Propanoate (ko00640) - succinate or acrylate route to propionate
#   - Butanoate (ko00650) - butyryl-CoA dehydrogenase consumes H2
#   - Pyruvate (ko00620) - central node for SCFA branching
#   - Carbon fixation prokaryotes (ko00720) - includes Wood-Ljungdahl acetogenesis
#   - Sulfur metabolism (ko00920) - dissimilatory sulfate reduction
#   - Nitrogen metabolism (ko00910) - dissimilatory nitrate reduction
#   - Glycolysis (ko00010) - acetate precursor flux
#
# Plus key KO-level H2-relevant enzymes from ANCOM-BC2 results.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

KEGG_FILE  <- "/mnt/d/Wisconsin_data1/merged_results/subset/kegg_pathway_matrix_subset.tsv"
ANCOM_FILE <- "results_ko/ancombc/ancombc_function_results.tsv"
OUT_DIR    <- "results_ko/h2_sink"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

meta <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_42","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45",
                  "RUM_KSU_46","RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","CON","SCFP","SCFP",
                  "CON","CON","SCFP")
)

# ---- A. Pathway level test ------------------------------------------------
pw <- fread(KEGG_FILE)
pw_cols <- setdiff(colnames(pw), "Sample")
ko_cols <- pw_cols[grepl("^ko", pw_cols)]
pw <- pw[, c("Sample", ko_cols), with = FALSE]
pw <- merge(meta, pw, by.x = "sample_id", by.y = "Sample")

# Drop zero-library samples
pw_mat <- as.matrix(pw[, ..ko_cols])
lib <- rowSums(pw_mat)
keep <- lib > 0
pw <- pw[keep,]
pw_mat <- pw_mat[keep,]
lib <- lib[keep]
pw_norm <- sweep(pw_mat, 1, lib, "/")
treat <- pw$treatment

focus <- c(
    "ko00680" = "Methane metabolism (reference, the H2 sink being suppressed)",
    "ko00640" = "Propanoate metabolism (succinate or acrylate route, classical H2 sink)",
    "ko00650" = "Butanoate metabolism (butyryl-CoA dehydrogenase, modest H2 sink)",
    "ko00620" = "Pyruvate metabolism (central SCFA branchpoint)",
    "ko00720" = "Carbon fixation in prokaryotes (includes Wood-Ljungdahl acetogenesis)",
    "ko00920" = "Sulfur metabolism (dissimilatory sulfate reduction H2 sink)",
    "ko00910" = "Nitrogen metabolism (dissimilatory nitrate reduction H2 sink)",
    "ko00010" = "Glycolysis (acetate precursor flux)",
    "ko00500" = "Starch and sucrose metabolism (substrate supply)",
    "ko00630" = "Glyoxylate and dicarboxylate metabolism (anaplerotic)"
)

pw_dt <- rbindlist(lapply(names(focus), function(p) {
    if (!(p %in% colnames(pw_norm))) return(NULL)
    v <- pw_norm[, p]
    con  <- v[treat == "CON"];  scfp <- v[treat == "SCFP"]
    wt <- suppressWarnings(wilcox.test(con, scfp))
    data.table(
        pathway = p, name = focus[[p]],
        CON_rel = signif(mean(con), 3),
        SCFP_rel = signif(mean(scfp), 3),
        log2_fc = signif(log2((mean(scfp) + 1e-9) / (mean(con) + 1e-9)), 3),
        direction = ifelse(mean(scfp) > mean(con), "UP in SCFP", "DOWN in SCFP"),
        W = unname(wt$statistic),
        p = wt$p.value
    )
}))
pw_dt[, q := p.adjust(p, method = "BH")]
setorder(pw_dt, p)
cat("=== Pathway level: alternative H2 sinks ===\n")
print(pw_dt)
fwrite(pw_dt, file.path(OUT_DIR, "h2_sink_pathway_compare.tsv"), sep = "\t")

# ---- B. KO level: specific H2-sink and SCFA enzymes ----------------------
anc <- fread(ANCOM_FILE)
h2_kos <- c(
    # Methanogenesis (suppressed - reference)
    "K00399" = "mcrA  methyl-CoM reductase alpha (methanogenesis)",
    "K00400" = "mcrC  methyl-CoM reductase C  (methanogenesis)",
    "K00401" = "mcrB  methyl-CoM reductase beta (methanogenesis)",
    "K00402" = "mcrG  methyl-CoM reductase gamma (methanogenesis)",
    "K14080" = "mtrA  N5-methyl-tetrahydromethanopterin (methanogenesis)",
    # Propionate via succinate pathway
    "K01847" = "MUT   methylmalonyl-CoA mutase (succinate to propionate)",
    "K00244" = "frdA  fumarate reductase A (succinate to propionate, H2 consuming)",
    "K00245" = "frdB  fumarate reductase B",
    "K00246" = "frdC  fumarate reductase C",
    "K00247" = "frdD  fumarate reductase D",
    "K01026" = "pct   propionate CoA-transferase",
    # Propionate via acrylate pathway (lactate to propionate)
    "K01595" = "ppc   PEP carboxylase",
    "K01596" = "pckA  PEP carboxykinase",
    # Butyrate
    "K00929" = "buk   butyrate kinase",
    "K00925" = "ackA  acetate kinase (acetate production)",
    "K00634" = "ptb   phosphate butyryltransferase",
    "K00248" = "bcd   butyryl-CoA dehydrogenase (H2 consuming)",
    # Lactate
    "K00016" = "ldh   L-lactate dehydrogenase",
    "K00101" = "lldD  L-lactate dehydrogenase quinone",
    # Wood-Ljungdahl reductive acetogenesis
    "K14138" = "acsB  acetyl-CoA synthase",
    "K00198" = "cooS  CO dehydrogenase",
    "K01938" = "fhs   formate-tetrahydrofolate ligase",
    "K00297" = "metF  methylenetetrahydrofolate reductase",
    "K01491" = "folD  methylenetetrahydrofolate dehydrogenase",
    # Sulfate reduction
    "K11180" = "dsrA  dissimilatory sulfite reductase alpha",
    "K11181" = "dsrB  dissimilatory sulfite reductase beta",
    "K00394" = "aprA  APS reductase alpha",
    "K00395" = "aprB  APS reductase beta",
    # Nitrate respiration
    "K00370" = "narG  nitrate reductase alpha",
    "K00374" = "narI  nitrate reductase gamma",
    "K03385" = "nrfA  cytochrome c nitrite reductase",
    "K00362" = "nirB  nitrite reductase NADH large",
    # Hydrogenases (H2 evolution vs uptake)
    "K00532" = "Fe hydrogenase (H2 evolving)",
    "K14107" = "ehbA  energy-converting hydrogenase B (Methanobrev)",
    "K14105" = "ehaN  energy-converting hydrogenase A (Methanobrev)"
)

ko_dt <- rbindlist(lapply(names(h2_kos), function(k) {
    r <- anc[feature == k]
    if (nrow(r) == 0) return(data.table(
        ko = k, gene_function = h2_kos[[k]],
        lfc = NA_real_, p = NA_real_, q = NA_real_, present = FALSE
    ))
    data.table(
        ko = k, gene_function = h2_kos[[k]],
        lfc = r$lfc, p = r$p, q = r$q, present = TRUE
    )
}))
setorder(ko_dt, p, na.last = TRUE)
cat("\n=== KO level: H2-sink and SCFA-producing enzymes (ANCOM-BC2) ===\n")
print(ko_dt)
fwrite(ko_dt, file.path(OUT_DIR, "h2_sink_KO_compare.tsv"), sep = "\t")

# ---- C. Pathway grouping summary ------------------------------------------
cat("\n=== Pathway grouping summary ===\n")
sink_groups <- list(
    "Methanogenesis (reference)"  = c("K00399","K00400","K00401","K00402","K14080"),
    "Propionate succinate route"  = c("K01847","K00244","K00245","K00246","K00247","K01026"),
    "Propionate acrylate route"   = c("K01595","K01596"),
    "Butyrate"                    = c("K00929","K00634","K00248"),
    "Acetate"                     = c("K00925"),
    "Lactate"                     = c("K00016","K00101"),
    "Wood-Ljungdahl acetogenesis" = c("K14138","K00198","K01938","K00297","K01491"),
    "Sulfate reduction"           = c("K11180","K11181","K00394","K00395"),
    "Nitrate respiration"         = c("K00370","K00374","K03385","K00362"),
    "Hydrogenases (methanogen)"   = c("K14105","K14107")
)
grp_dt <- rbindlist(lapply(names(sink_groups), function(g) {
    found <- intersect(sink_groups[[g]], anc$feature)
    if (length(found) == 0) return(data.table(
        group = g, n_KOs_found = 0L, mean_lfc = NA_real_,
        n_negative_LFC = NA_integer_, n_positive_LFC = NA_integer_
    ))
    sub <- anc[feature %in% found]
    data.table(
        group = g,
        n_KOs_found = nrow(sub),
        mean_lfc = round(mean(sub$lfc, na.rm = TRUE), 2),
        n_negative_LFC = sum(sub$lfc < 0, na.rm = TRUE),
        n_positive_LFC = sum(sub$lfc > 0, na.rm = TRUE)
    )
}))
print(grp_dt)
fwrite(grp_dt, file.path(OUT_DIR, "h2_sink_group_summary.tsv"), sep = "\t")

cat("\nOutputs in:", OUT_DIR, "\n")
