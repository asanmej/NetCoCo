# Author: Alejandro Santos Mejías
# Date last update: 2026-11-26
# Script name: edge_specific_threshold.R
# Aims: Umbral nulo específico por arista a partir del modelo nulo ya guardado (step_0).

# step_0 guarda, para cada permutación b, red_ebic_null[upper.tri(red_ebic_null)]
# uno detrás de otro. Por tanto el vector nulo se reorganiza en una matriz
# pares x permutaciones y cada fila es la distribución nula de un par concreto.
# El orden de los pares es el de upper.tri() con las columnas de df tras quitar
# las enfermedades excluidas, así que var_names debe ser colnames(df) en ese orden.
edge_specific_threshold <- function(null_values, var_names, prob = 0.95) {
  p <- length(var_names)
  m <- p * (p - 1) / 2
  B <- length(null_values) / m
  if (B != round(B)) {
    stop("La longitud del modelo nulo (", length(null_values),
         ") no es múltiplo de ", m, " pares: las variables no coinciden con el modelo nulo.")
  }
  null_mat <- matrix(abs(null_values), nrow = m, ncol = B)
  thr_vec <- apply(null_mat, 1, quantile, probs = prob, na.rm = TRUE, names = FALSE)

  thr <- matrix(0, p, p, dimnames = list(var_names, var_names))
  thr[upper.tri(thr)] <- thr_vec
  thr <- thr + t(thr)
  attr(thr, "n_perm") <- B
  thr
}
