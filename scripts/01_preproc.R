# =========================
# 01_preproc.R  (CFT project)
# Author: David Marcos Cuesta (based on Ben Falandays)
# Purpose: Load and clean Trial/ET/Click data, detect and time-lock beeps,
#          apply robust exclusions, resample to 120 Hz, and compute velocities.
# Outputs:
#   Data/ETdata_clean.Rdata(.csv)     # non-uniform time, quality flags kept
#   Data/Interpdata_clean.Rdata(.csv) # uniform 120 Hz, ready for windowing
# =========================

## ---- Parameters ----
options(digits.secs = 6)
FREQ_HZ      <- 120
DT_MS        <- 1000 / FREQ_HZ          # 8.333... ms
PAD_BLINK_N  <- 30                      # rows to mask before/after invalid samples
SHELF_BOUND  <- 1.5                     # plane bounds (meters)
MAX_TRIALS   <- 24
TRIAL_LOSS_T <- 0.10                    # >10% invalid → drop trial (unless beep)
SUBJ_LOSS_T  <- 0.20                    # >20% trials dropped → drop subject
SCORE_MIN    <- 6                       # drop trials with score < 6 (unless beep)
SUBS         <- seq(1, 31)              # participants to load

## ---- Libraries ----
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(dplyr)
  library(ggplot2); library(pracma); library(clock)
  library(slider)
})

## ---- Helpers ----
deg2rad <- function(x) x*pi/180
rad2deg <- function(x) x*180/pi

# ---- Paths ----
DATA_RAW       <- "data/raw"
DATA_PROCESSED <- "data/processed"
RESULTS_DIR    <- "results"
dir.create(DATA_PROCESSED, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR,    showWarnings = FALSE, recursive = TRUE)

## =========================
## 1) Load Trialdata
## =========================
trial_txt <- list.files(DATA_RAW, pattern = "^Trialdata_\\d+\\.txt$", full.names = TRUE)

if (length(trial_txt) > 0) {
  Trialdata <- purrr::map_dfr(trial_txt, \(f) data.table::fread(f, sep = ",", header = TRUE))
} else if (file.exists(file.path(DATA_RAW, "Trialdata.Rdata"))) {
  load(file.path(DATA_RAW, "Trialdata.Rdata"))  # crea objeto Trialdata
} else {
  stop("No Trialdata found in data/raw (expected Trialdata_*.txt or Trialdata.Rdata).")
}

Trialdata <- Trialdata %>%
  dplyr::filter(practiceNr == 3) %>%
  dplyr::mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>%
  dplyr::select(-practiceNr) %>%
  dplyr::select(order(colnames(.)))

# Parse timestamps
Trialdata <- Trialdata %>%
  dplyr::mutate(
    ymd = clock::date_parse(stringr::str_extract(trialStartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    trialStartTime  = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd),clock::get_month(ymd),clock::get_day(ymd),
                                                                 trialStartTime_h, trialStartTime_m, trialStartTime_s, trialStartTime_ms,
                                                                 subsecond_precision = "millisecond")),
    ymd = clock::date_parse(stringr::str_extract(targ1StartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    targ1StartTime  = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd),clock::get_month(ymd),clock::get_day(ymd),
                                                                 targ1StartTime_h, targ1StartTime_m, targ1StartTime_s, targ1StartTime_ms,
                                                                 subsecond_precision = "millisecond")),
    ymd = clock::date_parse(stringr::str_extract(trialFinishTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    trialFinishTime = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd),clock::get_month(ymd),clock::get_day(ymd),
                                                                 trialFinishTime_h, trialFinishTime_m, trialFinishTime_s, trialFinishTime_ms,
                                                                 subsecond_precision = "millisecond"))
  ) %>%
  dplyr::select(-tidyselect::matches("_(h|m|s|ms)$"), -ymd) %>%
  dplyr::select(order(colnames(.)))


## =========================
## 2) Load ETdata
## =========================
# Preferimos rawETdata.Rdata si existe (más rápido); si no, importamos los .txt
if (file.exists(file.path(DATA_RAW, "rawETdata.Rdata"))) {
  load(file.path(DATA_RAW, "rawETdata.Rdata"))  # objeto ETdata
} else {
  et_txt <- list.files(DATA_RAW, pattern = "^ETdata_\\d+\\.txt$", full.names = TRUE)
  if (length(et_txt) == 0) stop("No ETdata found in data/raw (expected ETdata_*.txt or rawETdata.Rdata).")
  ETdata <- purrr::map_dfr(et_txt, \(f) data.table::fread(f, sep = ",", header = TRUE)) %>%
    dplyr::filter(practiceNr == 3) %>%
    dplyr::mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>%
    dplyr::select(-practiceNr)
  # cache opcional para futuras corridas
  save(ETdata, file = file.path(DATA_RAW, "rawETdata.Rdata"))
}

# Parse SysTime, recomputa t y dt
ETdata <- ETdata %>%
  dplyr::mutate(
    ymd = clock::date_parse(stringr::str_extract(SysTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    SysTime = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd),clock::get_month(ymd),clock::get_day(ymd),
                                                         SysTime_h, SysTime_m, SysTime_s, SysTime_ms,
                                                         subsecond_precision = "millisecond"))
  ) %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::mutate(t = as.numeric(SysTime - SysTime[1]),
                dt = c(0, diff(t))) %>%
  dplyr::ungroup() %>%
  dplyr::select(-tidyselect::matches("SysTime_(h|m|s|ms)$"), -ymd, -time_stamp) %>%
  dplyr::select(order(colnames(.)))

# Glitch inicial más grande
rm_row <- which.max(ETdata$dt) - 1
if (length(rm_row) == 1 && rm_row > 0) {
  ETdata <- ETdata[-rm_row,] %>%
    dplyr::group_by(uniqueID) %>%
    dplyr::mutate(t = as.numeric(SysTime - SysTime[1]), dt = c(0, diff(t))) %>%
    dplyr::ungroup()
}
ETdata <- ETdata %>% dplyr::filter(!(uniqueID == "8_7" & t >= 33663))

# Merge Trialdata y anclas relativas al inicio de trial
ETdata <- dplyr::left_join(ETdata, Trialdata) %>% dplyr::select(order(colnames(.)))
ETdata <- ETdata %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::mutate(
    targ1StartTime  = as.numeric(targ1StartTime  - SysTime[1]),
    trialFinishTime = as.numeric(trialFinishTime - SysTime[1])
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(order(colnames(.)))


## =========================
## 3) Beep events (start vs detection)
## =========================
beeps <- ETdata %>%
  group_by(uniqueID) %>%
  mutate(
    beepStart  = as.integer(beepStatus == 1 & lag(beepStatus, default = 0) == 0),
    beepDetect = as.integer(beepStatus == 2 & lag(beepStatus, default = 0) == 1)
  ) %>%
  summarise(
    beepStartTime  = ifelse(any(beepStart==1),  t[which(beepStart==1)[1]],  NA_real_),
    beepDetectTime = ifelse(any(beepDetect==1), t[which(beepDetect==1)[1]], NA_real_)
  )

ETdata <- ETdata %>%
  left_join(beeps, by="uniqueID") %>%
  mutate(
    beepDetectionTrial = as.integer(!is.na(beepDetectTime)),
    t_beepStartTime     = t - beepStartTime,
    t_beepDetectionTime = t - beepDetectTime
  ) %>%
  select(order(colnames(.)))

## =========================
## 4) Plane projection is assumed pre-computed
## =========================
if (file.exists(file.path(DATA_RAW, "preProcessedETdata.Rdata"))) {
  load(file.path(DATA_RAW, "preProcessedETdata.Rdata"))  # ETdata con planeIntersect_*
}

## =========================
## 4) Plane projection (compute if missing, using your exact column names)
## =========================
need_proj <- !all(c("planeIntersect_x","planeIntersect_y") %in% names(ETdata))

if (need_proj) {
  message("→ Computing planeIntersect_x/y from raw gaze + HMD pose…")
  
  nm <- names(ETdata)
  
  # helper: pick first that exists from a list of candidates; else NA
  pick_from <- function(cands) {
    hit <- intersect(cands, nm)
    if (length(hit)) hit[1] else NA_character_
  }
  
  # === Mapear nombres EXACTOS de tu dataset, con fallback a regex ===
  col_goLx <- pick_from(c("gaze_origin_L.x(mm)", "gaze_origin_L.x.mm.", "gaze_origin_L_x_mm"))
  col_goLy <- pick_from(c("gaze_origin_L.y(mm)", "gaze_origin_L.y.mm.", "gaze_origin_L_y_mm"))
  col_goLz <- pick_from(c("gaze_origin_L.z(mm)", "gaze_origin_L.z.mm.", "gaze_origin_L_z_mm"))
  col_goRx <- pick_from(c("gaze_origin_R.x(mm)", "gaze_origin_R.x.mm.", "gaze_origin_R_x_mm"))
  col_goRy <- pick_from(c("gaze_origin_R.y(mm)", "gaze_origin_R.y.mm.", "gaze_origin_R_y_mm"))
  col_goRz <- pick_from(c("gaze_origin_R.z(mm)", "gaze_origin_R.z.mm.", "gaze_origin_R_z_mm"))
  
  col_gdLx <- pick_from(c("gaze_direct_L.x"))
  col_gdLy <- pick_from(c("gaze_direct_L.y"))
  col_gdLz <- pick_from(c("gaze_direct_L.z"))
  col_gdRx <- pick_from(c("gaze_direct_R.x"))
  col_gdRy <- pick_from(c("gaze_direct_R.y"))
  col_gdRz <- pick_from(c("gaze_direct_R.z"))
  
  col_rotX <- pick_from(c("camRot_x","camRot.x"))
  col_rotY <- pick_from(c("camRot_y","camRot.y"))
  col_rotZ <- pick_from(c("camRot_z","camRot.z"))
  col_posX <- pick_from(c("camPos_x","camPos.x"))
  col_posY <- pick_from(c("camPos_y","camPos.y"))
  col_posZ <- pick_from(c("camPos_z","camPos.z"))
  
  col_validL <- pick_from(c("eye_valid_L"))
  col_validR <- pick_from(c("eye_valid_R"))
  col_openL  <- pick_from(c("openness_L"))
  col_openR  <- pick_from(c("openness_R"))
  col_pupL   <- pick_from(c("pupil_diameter_L(mm)","pupil_diameter_L.mm."))
  col_pupR   <- pick_from(c("pupil_diameter_R(mm)","pupil_diameter_R.mm."))
  
  required <- c(col_goLx,col_goLy,col_goLz,col_goRx,col_goRy,col_goRz,
                col_gdLx,col_gdLy,col_gdLz,col_gdRx,col_gdRy,col_gdRz,
                col_rotX,col_rotY,col_rotZ,col_posX,col_posY,col_posZ)
  
  if (any(is.na(required))) {
    stop("Faltan columnas crudas necesarias (gaze_origin_*, gaze_direct_*, camRot_*, camPos_*). Revisa los nombres reales en names(ETdata).")
  }
  
  # ---- helpers matemáticos ----
  deg2rad <- function(x) x*pi/180
  wrap_pi <- function(a) ((a + pi) %% (2*pi)) - pi
  Rx <- function(th) matrix(c(1,0,0, 0,cos(th),-sin(th), 0,sin(th),cos(th)), 3,3, byrow=TRUE)
  Ry <- function(th) matrix(c(cos(th),0,sin(th), 0,1,0, -sin(th),0,cos(th)), 3,3, byrow=TRUE)
  Rz <- function(th) matrix(c(cos(th),-sin(th),0, sin(th),cos(th),0, 0,0,1), 3,3, byrow=TRUE)
  Rcompose <- function(x,y,z) Rz(z) %*% Rx(x) %*% Ry(y)  # z-x-y como en Ben
  
  # ray_plane_intersect <- function(p, v, n = c(0,0,-1), p0 = c(0,0,3.2)) {
  #   denom <- sum(v * n)
  #   if (abs(denom) < 1e-9) return(c(NA_real_, NA_real_, NA_real_))
  #   t <- sum((p0 - p) * n) / denom
  #   p + t * v
  # }
  ray_plane_intersect_safe <- function(p, v, n = c(0,0,-1), p0 = c(0,0,3.2)) {
    # Si hay NA/Inf en p o v, no intentamos nada
    if (any(!is.finite(p)) || any(!is.finite(v))) return(c(NA_real_, NA_real_, NA_real_))
    denom <- sum(v * n)
    # Si denom es NA/Inf o ~0, tampoco
    if (!is.finite(denom) || abs(denom) < 1e-9) return(c(NA_real_, NA_real_, NA_real_))
    t <- sum((p0 - p) * n) / denom
    if (!is.finite(t)) return(c(NA_real_, NA_real_, NA_real_))
    p + t * v
  }
  
  # ---- construir locales (mm→m, flip x), rotar/trasladar a mundo ----
  ETdata <- ETdata %>%
    mutate(
      goLx = - .data[[col_goLx]] * 0.001,  goLy =  .data[[col_goLy]] * 0.001,  goLz =  .data[[col_goLz]] * 0.001,
      goRx = - .data[[col_goRx]] * 0.001,  goRy =  .data[[col_goRy]] * 0.001,  goRz =  .data[[col_goRz]] * 0.001,
      gdLx = - .data[[col_gdLx]],         gdLy =  .data[[col_gdLy]],         gdLz =  .data[[col_gdLz]],
      gdRx = - .data[[col_gdRx]],         gdRy =  .data[[col_gdRy]],         gdRz =  .data[[col_gdRz]],
      rX   = wrap_pi(deg2rad(.data[[col_rotX]])),
      rY   = wrap_pi(deg2rad(.data[[col_rotY]])),
      rZ   = wrap_pi(deg2rad(.data[[col_rotZ]])),
      pX   = .data[[col_posX]], pY = .data[[col_posY]], pZ = .data[[col_posZ]]
    )
  
  # calidad/pesos por ojo (si faltan, no penalizamos)
  eye_valid_Lx <- if (!is.na(col_validL)) as.integer(ETdata[[col_validL]] == 31) else rep(1L, nrow(ETdata))
  eye_valid_Rx <- if (!is.na(col_validR)) as.integer(ETdata[[col_validR]] == 31) else rep(1L, nrow(ETdata))
  if (!is.na(col_openL)) eye_valid_Lx[ETdata[[col_openL]] < .2] <- 0L
  if (!is.na(col_openR)) eye_valid_Rx[ETdata[[col_openR]] < .2] <- 0L
  if (!is.na(col_pupL))  eye_valid_Lx[ETdata[[col_pupL]]  < 1]  <- 0L
  if (!is.na(col_pupR))  eye_valid_Rx[ETdata[[col_pupR]]  < 1]  <- 0L
  
  # normalizar dirs y pesos L/R
  dirL_mag <- sqrt(ETdata$gdLx^2 + ETdata$gdLy^2 + ETdata$gdLz^2)
  dirR_mag <- sqrt(ETdata$gdRx^2 + ETdata$gdRy^2 + ETdata$gdRz^2)
  ETdata$gdLx <- ETdata$gdLx / pmax(dirL_mag, 1e-9)
  ETdata$gdLy <- ETdata$gdLy / pmax(dirL_mag, 1e-9)
  ETdata$gdLz <- ETdata$gdLz / pmax(dirL_mag, 1e-9)
  ETdata$gdRx <- ETdata$gdRx / pmax(dirR_mag, 1e-9)
  ETdata$gdRy <- ETdata$gdRy / pmax(dirR_mag, 1e-9)
  ETdata$gdRz <- ETdata$gdRz / pmax(dirR_mag, 1e-9)
  
  wsum <- eye_valid_Lx + eye_valid_Rx
  Lw <- ifelse(wsum > 0, eye_valid_Lx / wsum, 0)
  Rw <- ifelse(wsum > 0, eye_valid_Rx / wsum, 0)
  
  # rotar/trasladar + proyectar (vectorizado por filas con vapply)
#   rot_and_project <- function(i) {
#     R <- Rcompose(ETdata$rX[i], ETdata$rY[i], ETdata$rZ[i])
#     
#     L_origin_w <- as.vector(R %*% c(ETdata$goLx[i], ETdata$goLy[i], ETdata$goLz[i])) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
#     R_origin_w <- as.vector(R %*% c(ETdata$goRx[i], ETdata$goRy[i], ETdata$goRz[i])) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
#     
#     L_dir_w <- as.vector(R %*% c(ETdata$gdLx[i], ETdata$gdLy[i], ETdata$gdLz[i]))
#     R_dir_w <- as.vector(R %*% c(ETdata$gdRx[i], ETdata$gdRy[i], ETdata$gdRz[i]))
#     
#     C_origin <- Lw[i]*L_origin_w + Rw[i]*R_origin_w
#     C_dir    <- Lw[i]*L_dir_w    + Rw[i]*R_dir_w
#     C_dir    <- C_dir / sqrt(sum(C_dir^2))
#     
#     ray_plane_intersect(C_origin, C_dir, n = c(0,0,-1), p0 = c(0,0,3.2))
#   }
#   
#   proj <- vapply(seq_len(nrow(ETdata)), rot_and_project, FUN.VALUE = numeric(3))
#   proj <- t(proj)
#   ETdata$planeIntersect_x <- proj[,1]
#   ETdata$planeIntersect_y <- proj[,2]
#   ETdata$planeIntersect_z <- proj[,3]
# }
  rot_and_project_safe <- function(i) {
    # cualquier NA en pose → NA
    if (!all(is.finite(ETdata$rX[i], ETdata$rY[i], ETdata$rZ[i],
                       ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i]))) {
      return(c(NA_real_, NA_real_, NA_real_))
    }
    R <- Rcompose(ETdata$rX[i], ETdata$rY[i], ETdata$rZ[i])
    
    # orígenes/direcciones ojo L/R (si falta algo → NA)
    Li <- c(ETdata$goLx[i], ETdata$goLy[i], ETdata$goLz[i])
    Ri <- c(ETdata$goRx[i], ETdata$goRy[i], ETdata$goRz[i])
    Ld <- c(ETdata$gdLx[i], ETdata$gdLy[i], ETdata$gdLz[i])
    Rd <- c(ETdata$gdRx[i], ETdata$gdRy[i], ETdata$gdRz[i])
    if (any(!is.finite(c(Li,Ri,Ld,Rd)))) return(c(NA_real_, NA_real_, NA_real_))
    
    L_origin_w <- as.vector(R %*% Li) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
    R_origin_w <- as.vector(R %*% Ri) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
    L_dir_w    <- as.vector(R %*% Ld)
    R_dir_w    <- as.vector(R %*% Rd)
    
    # pesos L/R ya calculados (Lw, Rw); si ambos 0, no hay dato útil
    if (!is.finite(Lw[i] + Rw[i]) || (Lw[i] + Rw[i]) <= 0) {
      return(c(NA_real_, NA_real_, NA_real_))
    }
    C_origin <- Lw[i]*L_origin_w + Rw[i]*R_origin_w
    C_dir    <- Lw[i]*L_dir_w    + Rw[i]*R_dir_w
    
    # normaliza dir; si norma ~0 o NA → NA
    nrm <- sqrt(sum(C_dir^2))
    if (!is.finite(nrm) || nrm < 1e-9) return(c(NA_real_, NA_real_, NA_real_))
    C_dir <- C_dir / nrm
    
    ray_plane_intersect_safe(C_origin, C_dir, n = c(0,0,-1), p0 = c(0,0,3.2))
  }
  
  # --- y la llamada que recorre filas, por esta ---
  proj <- vapply(seq_len(nrow(ETdata)), rot_and_project_safe, FUN.VALUE = numeric(3))
  proj <- t(proj)
  ETdata$planeIntersect_x <- proj[,1]
  ETdata$planeIntersect_y <- proj[,2]
  ETdata$planeIntersect_z <- proj[,3]
  # rot_and_project_safe <- function(i) {
  #   # cualquier NA en pose → NA
  #   if (!all(is.finite(ETdata$rX[i], ETdata$rY[i], ETdata$rZ[i],
  #                      ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i]))) {
  #     return(c(NA_real_, NA_real_, NA_real_))
  #   }
  #   R <- Rcompose(ETdata$rX[i], ETdata$rY[i], ETdata$rZ[i])
  #   
  #   # orígenes/direcciones ojo L/R (si falta algo → NA)
  #   Li <- c(ETdata$goLx[i], ETdata$goLy[i], ETdata$goLz[i])
  #   Ri <- c(ETdata$goRx[i], ETdata$goRy[i], ETdata$goRz[i])
  #   Ld <- c(ETdata$gdLx[i], ETdata$gdLy[i], ETdata$gdLz[i])
  #   Rd <- c(ETdata$gdRx[i], ETdata$gdRy[i], ETdata$gdRz[i])
  #   if (any(!is.finite(c(Li,Ri,Ld,Rd)))) return(c(NA_real_, NA_real_, NA_real_))
  #   
  #   L_origin_w <- as.vector(R %*% Li) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
  #   R_origin_w <- as.vector(R %*% Ri) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
  #   L_dir_w    <- as.vector(R %*% Ld)
  #   R_dir_w    <- as.vector(R %*% Rd)
  #   
  #   # pesos L/R ya calculados (Lw, Rw); si ambos 0, no hay dato útil
  #   if (!is.finite(Lw[i] + Rw[i]) || (Lw[i] + Rw[i]) <= 0) {
  #     return(c(NA_real_, NA_real_, NA_real_))
  #   }
  #   C_origin <- Lw[i]*L_origin_w + Rw[i]*R_origin_w
  #   C_dir    <- Lw[i]*L_dir_w    + Rw[i]*R_dir_w
  #   
  #   # normaliza dir; si norma ~0 o NA → NA
  #   nrm <- sqrt(sum(C_dir^2))
  #   if (!is.finite(nrm) || nrm < 1e-9) return(c(NA_real_, NA_real_, NA_real_))
  #   C_dir <- C_dir / nrm
  #   
  #   ray_plane_intersect_safe(C_origin, C_dir, n = c(0,0,-1), p0 = c(0,0,3.2))
  # }
  # 
  # # --- y la llamada que recorre filas, por esta ---
  # proj <- vapply(seq_len(nrow(ETdata)), rot_and_project_safe, FUN.VALUE = numeric(3))
  # proj <- t(proj)
  # ETdata$planeIntersect_x <- proj[,1]
  # ETdata$planeIntersect_y <- proj[,2]
  # ETdata$planeIntersect_z <- proj[,3]

if (!all(c("planeIntersect_x","planeIntersect_y") %in% names(ETdata))) {
  stop("No se pudieron crear planeIntersect_x/y. Revisa los nombres reales de columnas (gaze_origin_*, gaze_direct_*, camRot_*, camPos_*).")
}




## =========================
## 5) Invalid samples + padding (blink/edges) and exclusions
## =========================
ETdata <- ETdata %>%
  mutate(
    invalidData = ifelse(planeIntersect_x < -SHELF_BOUND | planeIntersect_x >  SHELF_BOUND |
                           planeIntersect_y < -SHELF_BOUND | planeIntersect_y >  SHELF_BOUND, 1, 0),
    invalidData = tidyr::replace_na(invalidData, 1)
  ) %>%
  group_by(uniqueID) %>%
  mutate(
    # row-based padding around invalid samples (stable, simple)
    invalidWindow = slider::slide_int(invalidData, ~ max(.x), .before = PAD_BLINK_N, .after = PAD_BLINK_N)
  ) %>%
  ungroup()

## Trial-level loss
trial_loss <- ETdata %>% group_by(uniqueID) %>% summarise(propInvalid = mean(invalidWindow==1))
exclude_trials <- trial_loss %>% filter(propInvalid > TRIAL_LOSS_T) %>% pull(uniqueID)
# keep beep trials even if over threshold
keep_beeps <- ETdata %>% filter(uniqueID %in% exclude_trials, beepDetectionTrial==1) %>% pull(uniqueID) %>% unique()
ETdata <- ETdata %>% filter(!(uniqueID %in% setdiff(exclude_trials, keep_beeps)))

## Subject-level loss
trial_counts <- ETdata %>%
  distinct(uniqueID, subjectNr, trialNr) %>%
  count(subjectNr, name="n_kept") %>%
  mutate(propInvalidTrials = 1 - n_kept/MAX_TRIALS)
exclude_subs <- trial_counts %>% filter(propInvalidTrials > SUBJ_LOSS_T) %>% pull(subjectNr)
ETdata <- ETdata %>% filter(!(subjectNr %in% exclude_subs))

## Low-score trials (except beep trials)
scores <- ETdata %>% group_by(uniqueID) %>% slice(1) %>% select(uniqueID, score, beepDetectionTrial)
low_score_trials <- scores %>% filter(score <= SCORE_MIN, beepDetectionTrial==0) %>% pull(uniqueID)
ETdata <- ETdata %>% filter(!(uniqueID %in% low_score_trials))

## =========================
## 6) Velocities (raw, non-uniform time)
## =========================
vel_raw <- ETdata %>%
  filter(invalidWindow==0) %>%
  group_by(uniqueID) %>%
  mutate(
    dt_ms  = c(0, diff(t)),
    xdiff  = c(0, diff(planeIntersect_x)),
    ydiff  = c(0, diff(planeIntersect_y)),
    eucDist= sqrt(xdiff^2 + ydiff^2),
    euc_vel= eucDist / (dt_ms/1000),           # meters/second
    theta2D= rad2deg(atan2(ydiff, xdiff)),
    ang_vel= theta2D / (dt_ms/1000)            # degrees/second
  ) %>%
  ungroup() %>%
  select(uniqueID, t, dt_ms, euc_vel, ang_vel)

ETdata_clean <- ETdata %>%
  select(-dt) %>%
  left_join(vel_raw, by=c("uniqueID","t")) %>%
  mutate(t_targ1StartTime = t - targ1StartTime) %>%
  select(order(colnames(.)))

## =========================
## 7) Interpolation to a regular 120 Hz grid
## =========================
interp_one <- function(df){
  xi <- seq(0, max(df$t, na.rm=TRUE), by = DT_MS)
  px <- pracma::interp1(df$t, df$planeIntersect_x, xi, method="linear")
  py <- pracma::interp1(df$t, df$planeIntersect_y, xi, method="linear")
  out <- tibble(t = xi, planeIntersect_x = px, planeIntersect_y = py)
  out$subjectNr <- df$subjectNr[1]; out$trialNr <- df$trialNr[1]; out$uniqueID <- df$uniqueID[1]
  out
}

Interpdata <- ETdata_clean %>%
  distinct(uniqueID, subjectNr, trialNr) %>%
  pull(uniqueID) %>%
  purrr::map_dfr(\(id){
    cur <- ETdata_clean %>% filter(uniqueID==id)
    suppressWarnings(interp_one(cur))
  })

## Add anchors and fixed dt; compute velocities on uniform grid
Interpdata <- Interpdata %>%
  left_join(
    ETdata_clean %>% group_by(uniqueID) %>%
      summarise(
        targ1StartTime   = first(targ1StartTime),
        trialFinishTime  = first(trialFinishTime),
        beepStartTime    = first(beepStartTime),
        beepDetectTime   = first(beepDetectTime),
        trialID          = first(trialID),
        target0          = first(target0),
        target1          = first(target1),
        totalMoves       = first(totalMoves),
        correctMoves     = first(correctMoves),
        completion       = first(completion),
        score            = first(score),
        beepDetectionTrial = first(beepDetectionTrial)
      ),
    by = "uniqueID"
  ) %>%
  group_by(uniqueID) %>%
  mutate(
    t_targ1StartTime     = t - targ1StartTime,
    t_beepStartTime      = t - beepStartTime,
    t_beepDetectionTime  = t - beepDetectTime,
    dt = DT_MS
  ) %>%
  mutate(
    xdiff  = c(0, diff(planeIntersect_x)),
    ydiff  = c(0, diff(planeIntersect_y)),
    eucDist= sqrt(xdiff^2 + ydiff^2),
    euc_vel= eucDist / (dt/1000),        # m/s
    theta2D= rad2deg(atan2(ydiff, xdiff)),
    ang_vel= theta2D / (dt/1000)         # deg/s
  ) %>%
  select(-xdiff, -ydiff, -eucDist, -theta2D) %>%
  ungroup() %>%
  select(order(colnames(.)))

## =========================
## 8) Save outputs + exclusion log
## =========================
save(ETdata_clean, file = file.path(DATA_PROCESSED, "ETdata_clean.Rdata"))
save(Interpdata,  file = file.path(DATA_PROCESSED, "Interpdata_clean.Rdata"))
write.csv(ETdata_clean, file = file.path(DATA_PROCESSED, "ETdata_clean.csv"), row.names = FALSE)
write.csv(Interpdata,  file = file.path(DATA_PROCESSED, "Interpdata_clean.csv"), row.names = FALSE)

exclusionTally <- list(
  trials_dropped_nonbeep = length(setdiff(exclude_trials, keep_beeps)),
  trials_kept_beep       = length(keep_beeps),
  subjects_dropped       = length(exclude_subs),
  trials_dropped_lowscore= length(low_score_trials)
)
writeLines(capture.output(str(exclusionTally)), con = file.path(RESULTS_DIR, "exclusion_tally.txt"))

message("✅ Preprocessing finished. Outputs in data/processed/")

