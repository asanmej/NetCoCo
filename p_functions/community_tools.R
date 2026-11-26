# Author: Alejandro Santos Mejías
# Date last update: 2026-11-26
# Script name: community_tools.R
# Aims: Ejecución repetida, estabilidad (AMI), matriz de co-asignación, partición
#       de consenso y selección de gamma, comunes a Leiden, Spinglass_pos y Spinglass_neg.

# Etiquetas canónicas (orden de aparición): particiones iguales -> vectores idénticos
relabel <- function(m) {
  out <- match(m, unique(m))
  names(out) <- names(m)
  out
}

safe_ami <- function(a, b) {
  a <- relabel(a); b <- relabel(b)
  if (identical(unname(a), unname(b))) return(1)   # evita NaN con dos particiones triviales
  aricode::AMI(a, b)
}

# Spinglass exige grafo conexo: se ejecuta por componente conexa.
# Componentes de 1 nodo -> singleton; de 2 nodos -> juntos si la arista es positiva.
spinglass_by_component <- function(g, wattr, gamma, gamma_minus, implementation, spins) {
  comp <- components(g)
  memb <- integer(vcount(g))
  next_id <- 0L
  for (cid in seq_len(comp$no)) {
    idx <- which(comp$membership == cid)
    if (length(idx) == 1) {
      next_id <- next_id + 1L
      memb[idx] <- next_id
      next
    }
    sg <- induced_subgraph(g, idx)
    w <- edge_attr(sg, wattr)
    if (length(idx) == 2) {
      if (sum(w) > 0) {
        next_id <- next_id + 1L
        memb[idx] <- next_id
      } else {
        memb[idx] <- next_id + 1:2
        next_id <- next_id + 2L
      }
      next
    }
    cl <- if (implementation == "neg") {
      cluster_spinglass(sg, weights = w, gamma = gamma, gamma.minus = gamma_minus,
                        implementation = "neg", spins = spins)
    } else {
      cluster_spinglass(sg, weights = w, gamma = gamma, spins = spins)
    }
    memb[idx] <- next_id + as.integer(membership(cl))
    next_id <- max(memb)
  }
  memb
}

# Una ejecución con semilla propia. Devuelve NA si el algoritmo falla.
run_community_once <- function(g, method, gamma, seed, spins = 25, gamma_minus = gamma) {
  set.seed(seed)
  memb <- tryCatch({
    if (method == "Leiden") {
      as.integer(membership(cluster_leiden(
        g, objective_function = "modularity", weights = E(g)$weight,
        resolution_parameter = gamma, n_iterations = -1)))
    } else if (method == "Spinglass_pos") {
      spinglass_by_component(g, "weight", gamma, gamma_minus, "orig", spins)
    } else if (method == "Spinglass_neg") {
      spinglass_by_component(g, "correlation", gamma, gamma_minus, "neg", spins)
    } else {
      stop("Método desconocido: ", method)
    }
  }, error = function(e) {
    message("Fallo ", method, " gamma=", gamma, " seed=", seed, ": ", conditionMessage(e))
    rep(NA_integer_, vcount(g))
  })
  if (anyNA(memb)) return(memb)
  relabel(memb)
}

# R ejecuciones con semillas seed+1..seed+R. Paraleliza si n_workers > 1 (plan() fuera).
run_replicates <- function(g, method, gamma, n_runs, seed, spins = 25,
                           gamma_minus = gamma, n_workers = 1) {
  seeds <- seed + seq_len(n_runs)
  f <- function(s) run_community_once(g, method, gamma, s, spins, gamma_minus)
  res <- if (n_workers > 1) {
    furrr::future_map(seeds, f, .options = furrr::furrr_options(seed = TRUE, packages = "igraph"))
  } else {
    lapply(seeds, f)
  }
  runs <- do.call(cbind, res)
  rownames(runs) <- V(g)$name
  colnames(runs) <- paste0("seed_", seeds)
  runs[, colSums(is.na(runs)) == 0, drop = FALSE]
}

pairwise_ami <- function(runs) {
  R <- ncol(runs)
  A <- matrix(NA_real_, R, R)
  for (i in seq_len(R - 1)) {
    for (j in (i + 1):R) {
      A[i, j] <- A[j, i] <- safe_ami(runs[, i], runs[, j])
    }
  }
  A
}

# D[i, j] = proporción de ejecuciones en que i y j caen en la misma comunidad
coassignment_matrix <- function(runs) {
  n <- nrow(runs)
  D <- matrix(0, n, n, dimnames = list(rownames(runs), rownames(runs)))
  for (r in seq_len(ncol(runs))) D <- D + outer(runs[, r], runs[, r], "==")
  D / ncol(runs)
}

# Consenso idéntico para los tres métodos (no favorece a ninguno):
# enlace medio sobre 1 - D, corte en co-asignación media >= tau
consensus_from_coassignment <- function(D, tau = 0.5) {
  hc <- hclust(as.dist(1 - D), method = "average")
  relabel(cutree(hc, h = 1 - tau))
}

# Confianza: co-asignación media con su comunidad; en singletons, 1 - máx. co-asignación
consensus_confidence <- function(D, cons) {
  sapply(seq_along(cons), function(i) {
    mates <- setdiff(which(cons == cons[i]), i)
    if (length(mates) == 0) 1 - max(D[i, -i]) else mean(D[i, mates])
  })
}

# Analiza un (método, gamma): ejecuciones, estabilidad, consenso y resumen
analyse_gamma <- function(g, method, gamma, n_runs, seed, tau = 0.5, spins = 25,
                          gamma_minus = gamma, n_workers = 1) {
  runs <- run_replicates(g, method, gamma, n_runs, seed, spins, gamma_minus, n_workers)
  if (ncol(runs) < 2) {
    warning("Menos de 2 ejecuciones válidas: ", method, " gamma=", gamma)
    return(NULL)
  }
  A <- pairwise_ami(runs)
  D <- coassignment_matrix(runs)
  cons <- consensus_from_coassignment(D, tau)
  conf <- consensus_confidence(D, cons)
  medoid <- runs[, which.max(rowMeans(A, na.rm = TRUE))]

  # Singletons triviales: comunidades de 1 enfermedad que SÍ tiene aristas.
  # Los nodos aislados son estructurales (umbral), no dependen de gamma.
  iso <- degree(g) == 0
  tab <- table(cons)
  singleton_ids <- as.integer(names(tab)[tab == 1])
  ami_vals <- A[upper.tri(A)]

  summary_row <- tibble(
    method = method,
    gamma = gamma,
    gamma_minus = if (method == "Spinglass_neg") gamma_minus else NA_real_,
    n_runs = n_runs,
    n_runs_ok = ncol(runs),
    stability_mean = mean(ami_vals),
    stability_sd = sd(ami_vals),
    k_runs_median = median(apply(runs, 2, function(x) length(unique(x)))),
    k_consensus = length(unique(cons)),
    n_trivial = sum(cons %in% singleton_ids & !iso),
    n_isolates = sum(iso),
    ami_consensus_medoid = safe_ami(cons, medoid)
  )
  list(runs = runs, D = D, consensus = cons, confidence = conf,
       medoid = medoid, summary = summary_row)
}

# Selección de gamma afinado por método sobre el barrido.
# Elegible: k >= min_k y como mucho max_trivial singletons triviales.
#  - gamma_max: gamma elegible con máxima estabilidad media (máximo puntual).
#  - gamma_plateau: centro de la meseta más larga de gammas elegibles consecutivos cuyo
#    consenso apenas cambia (AMI entre consensos vecinos >= plateau_ami).
#    Empates en longitud -> mayor estabilidad media en la meseta.
select_gamma <- function(sweep_summary, sweep_consensus, min_k = 2, max_trivial = 1,
                         plateau_ami = 0.9) {
  s <- sweep_summary %>% arrange(gamma)
  n <- nrow(s)
  keys <- sprintf("%.2f", s$gamma)
  s$eligible <- s$k_consensus >= min_k & s$n_trivial <= max_trivial
  s$ami_next <- c(vapply(seq_len(n - 1), function(i)
    safe_ami(sweep_consensus[[keys[i]]], sweep_consensus[[keys[i + 1]]]), numeric(1)), NA)

  out <- tibble(method = s$method[1], gamma_max = NA_real_, gamma_plateau = NA_real_,
                plateau_from = NA_real_, plateau_to = NA_real_, plateau_len = 0L,
                n_eligible = sum(s$eligible))
  s$plateau <- FALSE
  if (!any(s$eligible)) return(list(selection = out, sweep = s))

  el <- which(s$eligible)
  out$gamma_max <- s$gamma[el[which.max(s$stability_mean[el])]]

  run_id <- rep(NA_integer_, n)
  cur <- 0L
  for (i in el) {
    linked <- i > 1 && isTRUE(s$eligible[i - 1]) && isTRUE(s$ami_next[i - 1] >= plateau_ami)
    if (!linked) cur <- cur + 1L
    run_id[i] <- cur
  }
  plateaus <- split(el, run_id[el])
  len <- vapply(plateaus, length, integer(1))
  stab <- vapply(plateaus, function(r) mean(s$stability_mean[r]), numeric(1))
  best <- plateaus[[order(-len, -stab)[1]]]
  s$plateau[best] <- TRUE

  out$gamma_plateau <- s$gamma[best[ceiling(length(best) / 2)]]
  out$plateau_from <- min(s$gamma[best])
  out$plateau_to <- max(s$gamma[best])
  out$plateau_len <- length(best)
  list(selection = out, sweep = s)
}
