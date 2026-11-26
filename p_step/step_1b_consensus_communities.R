# Author: Alejandro Santos Mejías
# Date last update: 2026-11-26
# Script name: step_1b_consensus_communities.R
# Aims: Red con umbral nulo específico por arista y comunidades de consenso para
#       Leiden, Spinglass_pos y Spinglass_neg en tres análisis:
#       (1) principal con gamma = 1, (2) barrido de sensibilidad de gamma y
#       (3) gamma afinado por método (máximo puntual de estabilidad y meseta).

rm(list = ls())
gc()

source("98_packages.R")

# ============================================================================
# PARÁMETROS
# ============================================================================
seed <- 29101997
gamma_ebic <- 0.25
null_prob <- 0.95                          # cuantil del nulo de CADA arista
gamma_main <- 1
gamma_grid <- round(seq(0.5, 2.0, by = 0.1), 2)
n_runs_main <- 200                         # análisis principal y gamma afinado
n_runs_sweep <- 100                        # cada punto del barrido
tau_consensus <- 0.5                       # corte de co-asignación del consenso
min_k <- 2                                 # mínimo de comunidades
max_trivial <- 1                           # máximo de comunidades de 1 enfermedad (con aristas)
plateau_ami <- 0.9                         # AMI mínimo entre consensos vecinos en la meseta
spins <- 25
n_workers <- max(1, future::availableCores() - 1)
methods <- c("Leiden", "Spinglass_pos", "Spinglass_neg")

# ============================================================================
# RUTAS
# ============================================================================
null_model_files <- list.files(paste0(path_output, "null_model"), full.names = TRUE)
path_out <- paste0(path_output, "community_consensus/", format(Sys.Date(), "%Y%m%d"), "/")
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

file_list <- list.files(path = path_data,
                        pattern = "(2012|2017|2022)_(H|M)_(60-69|70-79|80-89|90)\\.csv$",
                        full.names = TRUE)

if (n_workers > 1) plan(multisession, workers = n_workers)

write_partition <- function(res, nodes, path_dir) {
  if (!dir.exists(path_dir)) dir.create(path_dir, recursive = TRUE)
  out <- nodes %>%
    mutate(community = res$consensus[id],
           confidence = res$confidence[match(id, names(res$consensus))],
           community_medoid = res$medoid[id])
  write_csv(out, file.path(path_dir, "nodos_consenso.csv"))
  write.csv(res$D, file.path(path_dir, "coasignacion.csv"), fileEncoding = "UTF-8")
  saveRDS(res$runs, file.path(path_dir, "ejecuciones.rds"))
}

all_network <- list()
all_sweep <- list()
all_selection <- list()
all_final <- list()

for (file_current in file_list) {

  name_group <- tools::file_path_sans_ext(basename(file_current))
  key <- regmatches(name_group, regexpr("(2012|2017|2022)_(H|M)_(60-69|70-79|80-89|90)", name_group))
  null_file <- grep(paste0("_", key, "\\.csv$"), null_model_files, value = TRUE)
  if (length(null_file) != 1) {
    message("Sin modelo nulo único para ", name_group, ": se omite")
    next
  }
  cat("\n==== ", name_group, " ====\n")
  path_group <- paste0(path_out, name_group, "/")
  if (!dir.exists(path_group)) dir.create(path_group, recursive = TRUE)

  # --------------------------------------------------------------------------
  # Datos y red (igual que step_1)
  # --------------------------------------------------------------------------
  df <- read_csv(file_current, show_col_types = FALSE)
  rm_disease <- c("CHROMOSOMAL_ABNORMALITIES", "INFLAMMATORY_BOWEL_DISEASES",
                  "MULTIPLE_SCLEROSIS", "SCHIZOPHRENIA_AND_DELUSIONAL_DISEASES")
  if (grepl("_M_", name_group)) rm_disease <- c(rm_disease, "PROSTATE_DISEASES")
  df <- df %>% select(all_of(setdiff(names(df), rm_disease)))

  matriz_r <- tetrachoric(df, smooth = TRUE)$rho
  red_ebic <- EBICglasso(matriz_r, n = nrow(df), gamma = gamma_ebic)
  colnames(red_ebic) <- rownames(red_ebic) <- colnames(df)

  # --------------------------------------------------------------------------
  # Umbral específico por arista (mismo orden de columnas que step_0)
  # --------------------------------------------------------------------------
  null_rho <- read.csv(null_file)$x
  thr <- edge_specific_threshold(null_rho, colnames(df), prob = null_prob)
  thr_pooled <- quantile(abs(null_rho), null_prob, na.rm = TRUE, names = FALSE)

  edges_all <- as.data.frame(as.table(red_ebic)) %>%
    rename(from = Var1, to = Var2, correlation = Freq) %>%
    mutate(from = as.character(from), to = as.character(to)) %>%
    filter(from < to) %>%
    mutate(weight = abs(correlation),
           threshold_edge = thr[cbind(from, to)],
           kept_edge = weight > threshold_edge,
           kept_pooled = weight > thr_pooled)
  write_csv(edges_all, paste0(path_group, "aristas_todas_", name_group, ".csv"))

  edges_filter <- edges_all %>%
    filter(kept_edge) %>%
    select(from, to, correlation, weight, threshold_edge)

  prevalence <- colMeans(df, na.rm = TRUE)
  nodes <- data.frame(id = names(prevalence), label = names(prevalence),
                      prevalencia_pct = prevalence * 100, size = 10 + prevalence * 100)

  g_pos <- graph_from_data_frame(edges_filter %>% filter(correlation > 0),
                                 vertices = nodes, directed = FALSE)
  g_signed <- graph_from_data_frame(edges_filter, vertices = nodes, directed = FALSE)
  graphs <- list(Leiden = g_pos, Spinglass_pos = g_pos, Spinglass_neg = g_signed)

  n_pairs <- nrow(edges_all)
  all_network[[name_group]] <- tibble(
    stratum = name_group, n = nrow(df), n_nodes = nrow(nodes), n_pairs = n_pairs,
    n_perm = attr(thr, "n_perm"), threshold_pooled = thr_pooled,
    threshold_edge_median = median(edges_all$threshold_edge),
    threshold_edge_max = max(edges_all$threshold_edge),
    edges_pooled = sum(edges_all$kept_pooled), edges_edge = nrow(edges_filter),
    edges_pos = sum(edges_filter$correlation > 0), edges_neg = sum(edges_filter$correlation < 0),
    density_edge = nrow(edges_filter) / n_pairs,
    isolates_pos = sum(degree(g_pos) == 0), isolates_signed = sum(degree(g_signed) == 0)
  )
  print(all_network[[name_group]])

  # --------------------------------------------------------------------------
  # (2) Barrido de sensibilidad de gamma
  # --------------------------------------------------------------------------
  sweep_partitions <- list()
  for (m in methods) {
    sweep_rows <- list()
    sweep_cons <- list()
    for (gm in gamma_grid) {
      cat(m, "gamma =", gm, "\n")
      res <- analyse_gamma(graphs[[m]], m, gm, n_runs_sweep, seed, tau_consensus,
                           spins, gamma_minus = gm, n_workers = n_workers)
      if (is.null(res)) next
      k <- sprintf("%.2f", gm)
      sweep_rows[[k]] <- res$summary
      sweep_cons[[k]] <- res$consensus
    }
    sel <- select_gamma(bind_rows(sweep_rows), sweep_cons, min_k, max_trivial, plateau_ami)
    all_sweep[[paste(name_group, m)]] <- sel$sweep %>% mutate(stratum = name_group, .before = 1)
    all_selection[[paste(name_group, m)]] <- sel$selection %>% mutate(stratum = name_group, .before = 1)
    sweep_partitions[[m]] <- sweep_cons
  }
  saveRDS(sweep_partitions, paste0(path_group, "barrido_consensos_", name_group, ".rds"))
  write_csv(bind_rows(all_sweep[paste(name_group, methods)]),
            paste0(path_group, "barrido_resumen_", name_group, ".csv"))

  # --------------------------------------------------------------------------
  # (1) Principal gamma = 1 y (3) gamma afinado (máximo y meseta)
  # --------------------------------------------------------------------------
  for (m in methods) {
    sel <- all_selection[[paste(name_group, m)]]
    analyses <- c(main_gamma1 = gamma_main, tuned_max = sel$gamma_max,
                  tuned_plateau = sel$gamma_plateau)
    cache <- list()
    for (a in names(analyses)) {
      gm <- unname(analyses[a])
      if (is.na(gm)) {
        message(m, " ", a, ": ningún gamma elegible")
        next
      }
      k <- sprintf("%.2f", gm)
      if (is.null(cache[[k]])) {
        cat(m, a, "gamma =", gm, "\n")
        cache[[k]] <- analyse_gamma(graphs[[m]], m, gm, n_runs_main, seed, tau_consensus,
                                    spins, gamma_minus = gm, n_workers = n_workers)
      }
      res <- cache[[k]]
      if (is.null(res)) next
      write_partition(res, nodes, paste0(path_group, a, "/", m, "/"))
      all_final[[paste(name_group, m, a)]] <- res$summary %>%
        mutate(stratum = name_group, analysis = a, .before = 1)
    }
  }
}

if (n_workers > 1) plan(sequential)

write_csv(bind_rows(all_network), paste0(path_out, "resumen_redes.csv"))
write_csv(bind_rows(all_sweep), paste0(path_out, "resumen_barrido.csv"))
write_csv(bind_rows(all_selection), paste0(path_out, "seleccion_gamma.csv"))
write_csv(bind_rows(all_final), paste0(path_out, "resumen_particiones.csv"))
cat("Resultados guardados en:", path_out, "\n")
