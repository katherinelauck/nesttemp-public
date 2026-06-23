## ================================================================
## HELPER FUNCTIONS
## ================================================================

scale_cols <- function(data, cols) {
  data |> mutate(across(all_of(cols), ~ scale(.x)[, 1], .names = "{.col}_scaled"))
}
# Filter dat to species, drop NA rows for temp_var and min_temp_var, scale
# cols_to_scale, and append a squared column for temp_var (add_sq = TRUE).
prep_model_data <- function(dat, species, temp_var, min_temp_var, cols_to_scale,
                            add_sq = TRUE) {
  d <- dat |>
    dplyr::filter(
      Species == species,
      !is.na(.data[[temp_var]]),
      !is.na(.data[[min_temp_var]])
    ) |>
    scale_cols(cols_to_scale)
  if (add_sq) {
    sq_col <- paste0(temp_var, "_scaled_sq")
    d[[sq_col]] <- d[[paste0(temp_var, "_scaled")]]^2
  }
  d
}

make_temp_trans <- function(data, col) {
  m <- mean(data[[col]], na.rm = TRUE)
  s <- sd(data[[col]], na.rm = TRUE)
  list(
    mean = m, sd = s,
    trans = scales::trans_new(col,
      transform = function(x) (x - m) / s,
      inverse   = function(x) x * s + m
    )
  )
}

anova_int_tab <- function(m_full, m_addmin, m_addmax, m_noint, digits = 4) {
  if (inherits(m_full, "glm") && !inherits(m_full, "lmerMod")) {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint, test = "Chisq") |>
        tibble() |>
        rename(Chisq = Deviance) |>
        mutate(AIC = c(AIC(m_full), AIC(m_single), AIC(m_noint)), .before = Df) |>
        mutate(Model = c("no interaction", "single interaction", "both interacting"), .before = Df)
    }
    p_col <- "Pr(>Chi)"
    rnd_cols <- "Chisq"
  } else {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint) |>
        tibble() |>
        mutate(Model = c("no interaction", "single interaction", "both interacting"), .before = npar)
    }
    p_col <- "Pr(>Chisq)"
    rnd_cols <- c("AIC", "Chisq")
  }
  bind_rows(mk_row(m_addmin), mk_row(m_addmax)) |>
    as_tibble() |>
    mutate(across(where(is.numeric), ~ round(.x, digits = digits)),
      P = .data[[p_col]]
    ) |>
    dplyr::select(Model, AIC, Chisq, P) |>
    mutate(
      max_or_min = c(rep("Max temp", times = 3), rep("Min temp", times = 3)),
      across(all_of(rnd_cols), ~ round(.x, digits = 2)),
      P = if_else(P < 0.001, "<0.001", as.character(P))
    ) |>
    group_by(max_or_min) |>
    gt()
}

format_emmeans_table <- function(model, spec = "habitat", digits = 2,
                                  p_digits = 3, use_regrid = TRUE) {
  em <- emmeans(model, spec)
  if (use_regrid) em <- regrid(em)
  em |> pairs() |> as_tibble() |>
    mutate(
      across(c(where(is.numeric), -any_of("p.value")), ~ round(.x, digits = digits)),
      across(p.value, ~ round(.x, digits = p_digits))
    ) |>
    gt()
}

format_emtrends_table <- function(model, temp_var, trend_label = "Max temp trend",
                                   use_pairwise = FALSE, digits = 3) {
  trend_col <- paste0(temp_var, ".trend")
  ratio_col <- if (use_pairwise) "z.ratio" else "t.ratio"
  spec <- if (use_pairwise) pairwise ~ habitat else ~ habitat
  result <- emtrends(model, specs = spec, var = temp_var) |> test()
  if (use_pairwise) result <- result |> pluck("emtrends")
  rename_vec <- c(
    "Habitat" = "habitat",
    setNames(trend_col, trend_label),
    "Df" = "df",
    "T-ratio" = ratio_col,
    "P" = "p.value"
  )
  result |>
    mutate(
      across(where(is.numeric), ~ round(.x, digits = digits)),
      df = round(df),
      p.value = if_else(p.value == 0.000, "<0.001", as.character(p.value))
    ) |>
    rename(all_of(rename_vec)) |>
    gt()
}
