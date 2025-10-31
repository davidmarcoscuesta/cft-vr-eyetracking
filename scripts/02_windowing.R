# =========================
# 02_windowing.R  (CFT project)
# Author: David M. Cuesta + ChatGPT
# Purpose: Create sliding/event-locked windows on Interpdata_clean and compute features.
# Inputs:  data/processed/Interpdata_clean.Rdata (.csv as fallback)
# Outputs: data/windowed/Interpdata_windowed_<params>.Rdata / .csv
# =========================

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readr)
  library(slider)
  library(stringr)
  library(glue)
})

## ---- Paths ----
DATA_PROCESSED <- "data/processed"
DATA_WINDOWED  <- "data/windowed"
dir.create(DATA_WINDOWED, showWarnings = FALSE, recursive = TRUE)

## ---- Parameters (edit here) ----
# Temporal grid comes from Interpdata_clean: t_ms (ms), dt_ms (ms, constante)
WIN_MS        <- 400L       # tamaño de ventana (ms)
STEP_MS       <- 40L        # paso (ms) para centros
MODE          <- "event"    # "sliding" | "event"
ALIGN         <- "targ1"    # "targ1" | "beepStart" | "beepDetect" (solo si MODE=="event")
PRE_MS        <- -2000L     # rango antes del evento (solo modo event)
POST_MS       <- 2000L      # rango después del evento (solo modo event)
MIN_PROP_OK   <- 0.60       # % mínimo de muestras válidas dentro de la ventana
ADD_GEOM      <- TRUE       # calcular distancia recorrida y desplazamientos medios
ADD_SPEED     <- TRUE       # calcular stats sobre euc_vel
ADD_ANG       <- TRUE       # calcular stats sobre ang_vel
VERBOSE       <- TRUE

## ---- Helpers ----
rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))

# Calcula features robustas para una serie (valores no-NA)
feat_stats <- function(x) {
  c(
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    iqr    = IQR(x, na.rm = TRUE),
    min    = suppressWarnings(min(x, na.rm = TRUE)),
    max    = suppressWarnings(max(x, na.rm = TRUE)),
    rms    = rms(x),
    n      = sum(is.finite(x))
  )
}

# Dado el df de un trial (uniqueID), construye centros de ventana
centers_for_mode <- function(df, mode, step_ms, pre_ms, post_ms, align) {
  tvec <- df$t_ms
  if (mode == "sliding") {
    # Deslizamos por todo el trial
    rng <- range(tvec, na.rm = TRUE)
    seq(from = ceiling(rng[1]), to = floor(rng[2]), by = step_ms)
  } else {
    # Event-locked: usamos columna relativa según ancla
    idx_col <- switch(
      align,
      "targ1"     = "t_targ1StartTime_ms",
      "beepStart" = "t_beepStartTime_ms",
      "beepDetect"= "t_beepDetectionTime_ms",
      stop("ALIGN debe ser 'targ1' | 'beepStart' | 'beepDetect'")
    )
    if (!(idx_col %in% names(df))) return(numeric(0))
    rel <- df[[idx_col]]
    # centros en coordenadas relativas [PRE_MS, POST_MS]
    seq(from = pre_ms, to = post_ms, by = step_ms)
  }
}

# Devuelve sub-dataframe dentro de [center - WIN_MS/2, center + WIN_MS/2]
slice_window <- function(df, center, win_ms, mode, align) {
  half <- win_ms / 2
  if (mode == "sliding") {
    win_start <- center - half
    win_end   <- center + half
    df %>% dplyr::filter(t_ms >= win_start, t_ms <= win_end) %>%
      mutate(win_start_ms = win_start, win_end_ms = win_end, center_ms = center)
  } else {
    idx_col <- switch(
      align,
      "targ1"     = "t_targ1StartTime_ms",
      "beepStart" = "t_beepStartTime_ms",
      "beepDetect"= "t_beepDetectionTime_ms",
      stop("ALIGN debe ser 'targ1' | 'beepStart' | 'beepDetect'")
    )
    win_start <- center - half
    win_end   <- center + half
    df %>% dplyr::filter(.data[[idx_col]] >= win_start, .data[[idx_col]] <= win_end) %>%
      mutate(win_start_ms = win_start, win_end_ms = win_end, center_ms = center)
  }
}

# Métricas geométricas (opcional) usando el grid regular
geom_metrics <- function(win_df) {
  # desplazamientos (x,y) y distancia recorrida aproximada en la ventana
  # Nota: dt_ms es constante en Interpdata_clean
  out <- c()
  if (all(c("planeIntersect_x","planeIntersect_y","dt_ms") %in% names(win_df))) {
    x <- win_df$planeIntersect_x
    y <- win_df$planeIntersect_y
    dx <- c(0, diff(x))
    dy <- c(0, diff(y))
    step_dist <- sqrt(dx^2 + dy^2)       # en metros por paso
    # tiempo por paso
    dt <- if ("dt_ms" %in% names(win_df)) (win_df$dt_ms[1]/1000) else NA_real_
    total_dist <- sum(step_dist, na.rm = TRUE)
    mean_speed_geom <- if (is.finite(dt)) mean(step_dist / dt, na.rm = TRUE) else NA_real_
    out <- c(
      disp_x_mean = mean(dx, na.rm = TRUE),
      disp_y_mean = mean(dy, na.rm = TRUE),
      path_len    = total_dist,
      mean_speed_geom = mean_speed_geom
    )
  }
  out
}

# Resumen por ventana
summarise_window <- function(win_df, add_geom = TRUE, add_speed = TRUE, add_ang = TRUE) {
  n_all <- nrow(win_df)
  # Si no hay dt_ms en Interpdata, lo marcamos como NA (pero en nuestro pipeline sí existe)
  dt_s  <- if ("dt_ms" %in% names(win_df) && n_all > 0) win_df$dt_ms[1] / 1000 else NA_real_
  
  feats <- c()
  
  if (add_speed && "euc_vel" %in% names(win_df)) {
    sp <- win_df$euc_vel
    feats <- c(feats,
               setNames(feat_stats(sp), paste0("euc_vel_", names(feat_stats(sp)))))
  }
  
  if (add_ang && "ang_vel" %in% names(win_df)) {
    av <- win_df$ang_vel
    feats <- c(feats,
               setNames(feat_stats(av), paste0("ang_vel_", names(feat_stats(av)))))
  }
  
  if (add_geom) {
    feats <- c(feats, geom_metrics(win_df))
  }
  
  # conteos
  feats <- c(
    feats,
    n_samples = n_all
  )
  
  # Devuelve named numeric vector
  feats
}

# Control de calidad: suficiente proporción de muestras no-NA
window_ok <- function(win_df, win_ms, min_prop_ok) {
  if (nrow(win_df) == 0) return(FALSE)
  # esperado (aprox) según dt_ms
  if (!"dt_ms" %in% names(win_df)) return(TRUE)  # si no está, no filtramos
  expected <- ceiling(win_ms / win_df$dt_ms[1])
  valid_n  <- sum(is.finite(win_df$euc_vel) | is.finite(win_df$ang_vel) |
                    (is.finite(win_df$planeIntersect_x) & is.finite(win_df$planeIntersect_y)))
  (valid_n / expected) >= min_prop_ok
}

## ---- Load Interpdata ----
RDATA <- file.path(DATA_PROCESSED, "Interpdata_clean.Rdata")
CSV   <- file.path(DATA_PROCESSED, "Interpdata_clean.csv")

if (file.exists(RDATA)) {
  load(RDATA)  # expects object Interpdata
} else if (file.exists(CSV)) {
  Interpdata <- readr::read_csv(CSV, show_col_types = FALSE)
} else {
  stop("No Interpdata_clean found in data/processed/")
}

need_cols <- c("uniqueID","subjectNr","trialNr","t_ms","dt_ms",
               "planeIntersect_x","planeIntersect_y","euc_vel","ang_vel",
               "t_targ1StartTime_ms","t_beepStartTime_ms","t_beepDetectionTime_ms")
missing_cols <- setdiff(need_cols, names(Interpdata))
if (length(missing_cols)) {
  stop("Missing required columns in Interpdata_clean: ", paste(missing_cols, collapse = ", "))
}

## ---- Build windows per trial ----
if (VERBOSE) {
  message(glue("→ Windowing mode: {MODE} | align={ALIGN} | win={WIN_MS}ms | step={STEP_MS}ms | pre={PRE_MS} | post={POST_MS}"))
}

windows_df <- Interpdata %>%
  group_by(uniqueID, subjectNr, trialNr) %>%
  group_modify(function(df, key) {
    df <- arrange(df, t_ms)
    
    centers <- centers_for_mode(df, MODE, STEP_MS, PRE_MS, POST_MS, ALIGN)
    if (!length(centers)) return(tibble())
    
    # Para event mode, los centros están en tiempo relativo: generamos ventanas en coords relativas
    res_list <- vector("list", length(centers))
    for (i in seq_along(centers)) {
      wdf <- slice_window(df, centers[i], WIN_MS, MODE, ALIGN)
      ok  <- window_ok(wdf, WIN_MS, MIN_PROP_OK)
      if (!ok) next
      
      feats <- summarise_window(wdf, add_geom = ADD_GEOM, add_speed = ADD_SPEED, add_ang = ADD_ANG)
      res_list[[i]] <- tibble(
        center_ms   = wdf$center_ms[1],
        win_start_ms= wdf$win_start_ms[1],
        win_end_ms  = wdf$win_end_ms[1],
        mode        = MODE,
        align       = if (MODE=="event") ALIGN else "none",
        !!!as.list(feats)
      )
    }
    bind_rows(res_list)
  }) %>%
  ungroup()

if (nrow(windows_df) == 0) {
  warning("No windows produced (check parameters or MIN_PROP_OK too strict).")
}

## ---- Optional: pick closest-to-zero windows around targ1 (one negative & one positive) ----
closest_zero <- function(df, center_col = "center_ms") {
  # negativo más cercano a 0 y positivo más cercano a 0
  neg <- df %>% filter(.data[[center_col]] <= 0) %>% slice_max(.data[[center_col]], n = 1, with_ties = FALSE)
  pos <- df %>% filter(.data[[center_col]] >= 0) %>% slice_min(.data[[center_col]], n = 1, with_ties = FALSE)
  bind_rows(neg, pos)
}

closest_tbl <- tibble()
if (MODE == "event" && ALIGN == "targ1" && nrow(windows_df)) {
  closest_tbl <- windows_df %>%
    group_by(uniqueID) %>%
    group_modify(~ closest_zero(.x, "center_ms")) %>%
    ungroup() %>%
    mutate(tag = "closest_to_0")
}

## ---- Save ----
param_tag <- if (MODE == "sliding") {
  glue("win{WIN_MS}_step{STEP_MS}_mode-sliding")
} else {
  glue("win{WIN_MS}_step{STEP_MS}_mode-event_align-{ALIGN}_pre{PRE_MS}_post{POST_MS}")
}

rds_path  <- file.path(DATA_WINDOWED, glue("Interpdata_windowed_{param_tag}.Rdata"))
csv_path  <- file.path(DATA_WINDOWED, glue("Interpdata_windowed_{param_tag}.csv"))
rds_close <- file.path(DATA_WINDOWED, glue("Interpdata_windowed_closest0_{param_tag}.Rdata"))
csv_close <- file.path(DATA_WINDOWED, glue("Interpdata_windowed_closest0_{param_tag}.csv"))

save(windows_df, file = rds_path)
readr::write_csv(windows_df, csv_path)

if (nrow(closest_tbl)) {
  save(closest_tbl, file = rds_close)
  readr::write_csv(closest_tbl, csv_close)
}

if (VERBOSE) {
  message("✅ Windowing done.")
  message("→ Saved: ", rds_path)
  message("→ Saved: ", csv_path)
  if (nrow(closest_tbl)) {
    message("→ Saved closest-to-zero around targ1: ", rds_close)
    message("→ Saved closest-to-zero around targ1: ", csv_close)
  }
}