# simular_datos.R ---------------------------------------------------------
# Genera un dataset FICTICIO con la misma estructura que el del estudio.
# Los datos reales no se publican: son de pacientes.
#
# Qué reproduce (de forma aproximada):
#   - n = 1685 pacientes; solo 582 con edad registrada
#     (0-5: 44; 6-11: 73; 12-17: 70; 18+: 395).
#   - La prevalencia de sensibilización (IgE específica >= 0,35 kUA/L) a los
#     nueve componentes de Der p en cada grupo de edad.
#   - La co-sensibilización: Der p 1, 7 y 20 van muy juntos; Der p 2, 5 y 10
#     van casi por libre.
#   - La relación entre IgE total y número de componentes positivos.
#   - Algún valor ausente en los componentes, como en un informe real.
# Las cifras no serán idénticas a las publicadas.
#
# Uso: source("R/simular_datos.R") crea data/der_p_simulado.csv
# -------------------------------------------------------------------------

set.seed(2025)

allergens <- c("der_p_1", "der_p_2", "der_p_5", "der_p_7", "der_p_10",
               "der_p_11", "der_p_20", "der_p_21", "der_p_23")

# Tamaño de cada grupo de edad ("unknown" = sin fecha de nacimiento)
group_n <- c("0-5" = 44, "6-11" = 73, "12-17" = 70, "18+" = 395, "unknown" = 1103)

# Prevalencia objetivo (%) por grupo de edad (filas) y componente (columnas)
prev <- rbind(
  "0-5"     = c(38.6,  6.8, 0.0, 25.0, 2.3,  9.1, 22.7,  6.8,  4.5),
  "6-11"    = c(64.8, 23.9, 1.4, 53.3, 5.6, 30.4, 52.9, 42.5, 24.3),
  "12-17"   = c(70.0, 10.0, 0.0, 71.4, 7.1, 34.8, 70.0, 40.0, 39.1),
  "18+"     = c(34.9,  4.6, 0.3, 42.6, 5.1, 21.5, 39.7, 25.8, 20.8),
  "unknown" = c(48.2,  7.1, 0.2, 51.2, 7.9, 26.8, 48.0, 30.6, 27.3)
) / 100
colnames(prev) <- allergens

# Cuánto depende cada componente de la "atopia" común del paciente.
# Alto = muy co-sensibilizado con los demás; bajo = va por libre.
loading <- c(der_p_1 = 2.4, der_p_2 = 0.3, der_p_5 = 0.2, der_p_7 = 2.4,
             der_p_10 = 0.4, der_p_11 = 1.5, der_p_20 = 2.4, der_p_21 = 1.6,
             der_p_23 = 1.5)

# Desplazamiento de la IgE total por edad (log), con pico en adolescentes
tige_shift <- c("0-5" = -0.8, "6-11" = 0.2, "12-17" = 0.6, "18+" = 0, "unknown" = 0)

age_ranges <- list("0-5" = c(0, 5), "6-11" = c(6, 11), "12-17" = c(12, 17),
                   "18+" = c(18, 90))

sim_group <- function(g, n) {
  z <- rnorm(n)  # atopia latente de cada paciente
  out <- data.frame(age_group = rep(g, n))
  for (a in allergens) {
    p <- min(max(prev[g, a], 1e-4), 1 - 1e-4)
    lam <- loading[[a]]
    # Umbral de la variable latente que da la prevalencia buscada
    b <- qnorm(p) * sqrt(lam^2 + 1)
    positive <- lam * z + b + rnorm(n) > 0
    level_pos <- 0.35 + exp(rnorm(n, log(1.5) + 0.5 * z, 0.7))
    level_neg <- ifelse(runif(n) < 0.7, 0.10, runif(n, 0.10, 0.34))  # "< 0,10" se guarda como 0,10
    out[[a]] <- round(ifelse(positive, level_pos, level_neg), 2)
  }
  n_pos <- rowSums(out[allergens] >= 0.35)
  out$total_ige <- round(exp(log(40) + 0.20 * n_pos + 0.25 * z + tige_shift[[g]] + rnorm(n, 0, 1.1)), 1)
  out$age <- if (g == "unknown") NA_integer_ else
    sample(age_ranges[[g]][1]:age_ranges[[g]][2], n, replace = TRUE)
  out
}

der_p <- do.call(rbind, Map(sim_group, names(group_n), group_n))

# Algún valor ausente en los componentes (~3 %), como en un informe real
for (a in allergens) der_p[[a]][runif(nrow(der_p)) < 0.03] <- NA

der_p <- der_p[sample(nrow(der_p)), ]
der_p$id <- sprintf("P%04d", seq_len(nrow(der_p)))
der_p <- der_p[, c("id", "age", "total_ige", allergens)]

dir.create("data", showWarnings = FALSE)
write.csv(der_p, "data/der_p_simulado.csv", row.names = FALSE, na = "")
