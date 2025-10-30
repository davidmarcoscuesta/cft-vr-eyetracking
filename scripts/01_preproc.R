# =========================
# 01_preproc.R  (CFT project)
# Author: David Cuesta (based on Ben Falandays) + hardening
# Purpose: Load and clean Trial / ET data, detect and time-lock beeps,
#          apply exclusions, project gaze to shelf plane, resample to 120 Hz,
#          compute velocities, and write processed outputs.
# Outputs:
#   data/processed/ETdata_clean.Rdata(.csv)     # non-uniform time, quality flags kept
#   data/processed/Interpdata_clean.Rdata(.csv) # uniform 120 Hz, ready for windowing
#   results/exclusion_tally.txt
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
SUBS         <- seq(1, 31)              # fallback, but we autodetect files in data/raw

## ---- Libraries ----
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(pracma)
  library(clock)
  library(slider)
  library(stringr)
})

## ---- Helpers ----
deg2rad <- function(x) x*pi/180
rad2deg <- function(x) x*180/pi

## ---- Paths ----
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
} else if (file.exists(file.path(DATA_RAW, "Trialdata.csv"))) {
  Trialdata <- readr::read_csv(file.path(DATA_RAW, "Trialdata.csv"), show_col_types = FALSE)
} else if (file.exists(file.path(DATA_RAW, "Trialdata.Rdata"))) {
  load(file.path(DATA_RAW, "Trialdata.Rdata"))  # expects object Trialdata
} else {
  stop("No Trialdata found in data/raw (expected Trialdata_*.txt, Trialdata.csv or Trialdata.Rdata).")
}

Trialdata <- Trialdata %>%
  dplyr::filter(practiceNr == 3) %>%
  dplyr::mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>%
  dplyr::select(-practiceNr)

# Parse timestamps to clock::naive_time
Trialdata <- Trialdata %>%
  dplyr::mutate(
    ymd = clock::date_parse(stringr::str_extract(trialStartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    trialStartTime  = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd), clock::get_month(ymd), clock::get_day(ymd),
                                                                 trialStartTime_h, trialStartTime_m, trialStartTime_s, trialStartTime_ms,
                                                                 subsecond_precision = "millisecond")),
    ymd = clock::date_parse(stringr::str_extract(targ1StartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    targ1StartTime  = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd), clock::get_month(ymd), clock::get_day(ymd),
                                                                 targ1StartTime_h, targ1StartTime_m, targ1StartTime_s, targ1StartTime_ms,
                                                                 subsecond_precision = "millisecond")),
    ymd = clock::date_parse(stringr::str_extract(trialFinishTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    trialFinishTime = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd), clock::get_month(ymd), clock::get_day(ymd),
                                                                 trialFinishTime_h, trialFinishTime_m, trialFinishTime_s, trialFinishTime_ms,
                                                                 subsecond_precision = "millisecond"))
  ) %>%
  dplyr::select(-matches("_(h|m|s|ms)$"), -ymd) %>%
  dplyr::select(order(colnames(.)))

## =========================
## 2) Load ETdata
## =========================
# Prefer a raw Rdata dump if present (fast), else stitch ETdata_*.txt
if (file.exists(file.path(DATA_RAW, "rawETdata.Rdata"))) {
  load(file.path(DATA_RAW, "rawETdata.Rdata"))  # should create ETdata
  if (!exists("ETdata")) stop("rawETdata.Rdata loaded but ETdata object not found.")
} else {
  et_txt <- list.files(DATA_RAW, pattern = "^ETdata_\\d+\\.txt$", full.names = TRUE)
  if (length(et_txt) == 0) stop("No ETdata found in data/raw (expected ETdata_*.txt or rawETdata.Rdata).")
  ETdata <- purrr::map_dfr(et_txt, \(f) data.table::fread(f, sep = ",", header = TRUE))
}

## Parse SysTime, recompute t and dt
## Parse SysTime, recompute t_ms y dt_ms (ambos en milisegundos)
ETdata <- ETdata %>%
  dplyr::mutate(
    ymd = clock::date_parse(stringr::str_extract(SysTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    SysTime = clock::as_naive_time(clock::year_month_day(clock::get_year(ymd), clock::get_month(ymd), clock::get_day(ymd),
                                                         SysTime_h, SysTime_m, SysTime_s, SysTime_ms,
                                                         subsecond_precision = "millisecond"))
  ) %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::mutate(
    # t_ms: milisegundos desde el inicio del trial
    t_ms  = as.numeric((SysTime - SysTime[1]) / clock::duration_milliseconds(1)),
    dt_ms = c(0, diff(t_ms))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-matches("SysTime_(h|m|s|ms)$"), -ymd, -time_stamp) %>%
  dplyr::select(order(colnames(.)))

## Remove the single largest dt glitch and re-compute; known tail glitch in "8_7"
rm_row <- which.max(ETdata$dt) - 1
if (length(rm_row) == 1 && rm_row > 0) {
  ETdata <- ETdata[-rm_row,] %>%
    dplyr::group_by(uniqueID) %>%
    dplyr::mutate(t = as.numeric(SysTime - SysTime[1]), dt = c(0, diff(t))) %>%
    dplyr::ungroup()
}
ETdata <- ETdata %>% dplyr::filter(!(uniqueID == "8_7" & t >= 33663))

## Merge Trialdata
ETdata <- dplyr::left_join(ETdata, Trialdata, by = c("subjectNr","trialNr","uniqueID")) %>%
  dplyr::select(order(colnames(.)))

## Make trial-relative anchors (ms from trial start)
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
  dplyr::group_by(uniqueID) %>%
  dplyr::mutate(
    beepStart  = as.integer(beepStatus == 1 & dplyr::lag(beepStatus, default = 0) == 0),
    beepDetect = as.integer(beepStatus == 2 & dplyr::lag(beepStatus, default = 0) == 1)
  ) %>%
  dplyr::summarise(
    beepStartTime  = ifelse(any(beepStart==1),  t[which(beepStart==1)[1]],  NA_real_),
    beepDetectTime = ifelse(any(beepDetect==1), t[which(beepDetect==1)[1]], NA_real_)
  )

ETdata <- ETdata %>%
  dplyr::left_join(beeps, by = "uniqueID") %>%
  dplyr::mutate(
    beepDetectionTrial = as.integer(!is.na(beepDetectTime)),
    t_beepStartTime     = t - beepStartTime,
    t_beepDetectionTime = t - beepDetectTime
  ) %>%
  dplyr::select(order(colnames(.)))

## =========================
## 4) Plane projection (compute if missing, using your exact column names)
## =========================
need_proj <- !all(c("planeIntersect_x","planeIntersect_y") %in% names(ETdata))

if (need_proj) {
  message("→ Computing planeIntersect_x/y from raw gaze + HMD pose…")
  
  nm <- names(ETdata)
  pick_from <- function(cands) { hit <- intersect(cands, nm); if (length(hit)) hit[1] else NA_character_ }
  
  # Exact names in your files (with fallback variants)
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
    stop("Missing raw columns for projection (gaze_origin_*, gaze_direct_*, camRot_*, camPos_*).")
  }
  
  # math helpers
  wrap_pi <- function(a) ((a + pi) %% (2*pi)) - pi
  Rx <- function(th) matrix(c(1,0,0, 0,cos(th),-sin(th), 0,sin(th),cos(th)), 3,3, byrow=TRUE)
  Ry <- function(th) matrix(c(cos(th),0,sin(th), 0,1,0, -sin(th),0,cos(th)), 3,3, byrow=TRUE)
  Rz <- function(th) matrix(c(cos(th),-sin(th),0, sin(th),cos(th),0, 0,0,1), 3,3, byrow=TRUE)
  Rcompose <- function(x,y,z) Rz(z) %*% Rx(x) %*% Ry(y)  # z-x-y like Ben
  
  ray_plane_intersect_safe <- function(p, v, n = c(0,0,-1), p0 = c(0,0,3.2)) {
    if (any(!is.finite(p)) || any(!is.finite(v))) return(c(NA_real_, NA_real_, NA_real_))
    denom <- sum(v * n)
    if (!is.finite(denom) || abs(denom) < 1e-9) return(c(NA_real_, NA_real_, NA_real_))
    t <- sum((p0 - p) * n) / denom
    if (!is.finite(t)) return(c(NA_real_, NA_real_, NA_real_))
    p + t * v
  }
  
  ETdata <- ETdata %>%
    dplyr::mutate(
      goLx = - .data[[col_goLx]] * 0.001,  goLy =  .data[[col_goLy]] * 0.001,  goLz =  .data[[col_goLz]] * 0.001,
      goRx = - .data[[col_goRx]] * 0.001,  goRy =  .data[[col_goRy]] * 0.001,  goRz =  .data[[col_goRz]] * 0.001,
      gdLx = - .data[[col_gdLx]],         gdLy =  .data[[col_gdLy]],         gdLz =  .data[[col_gdLz]],
      gdRx = - .data[[col_gdRx]],         gdRy =  .data[[col_gdRy]],         gdRz =  .data[[col_gdRz]],
      rX   = wrap_pi(deg2rad(.data[[col_rotX]])),
      rY   = wrap_pi(deg2rad(.data[[col_rotY]])),
      rZ   = wrap_pi(deg2rad(.data[[col_rotZ]])),
      pX   = .data[[col_posX]], pY = .data[[col_posY]], pZ = .data[[col_posZ]]
    )
  
  # eye validity → weights
  eye_valid_Lx <- if (!is.na(col_validL)) as.integer(ETdata[[col_validL]] == 31) else rep(1L, nrow(ETdata))
  eye_valid_Rx <- if (!is.na(col_validR)) as.integer(ETdata[[col_validR]] == 31) else rep(1L, nrow(ETdata))
  if (!is.na(col_openL)) eye_valid_Lx[ETdata[[col_openL]] < .2] <- 0L
  if (!is.na(col_openR)) eye_valid_Rx[ETdata[[col_openR]] < .2] <- 0L
  if (!is.na(col_pupL))  eye_valid_Lx[ETdata[[col_pupL]]  < 1]  <- 0L
  if (!is.na(col_pupR))  eye_valid_Rx[ETdata[[col_pupR]]  < 1]  <- 0L
  
  # normalize directions and compute weights
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
  
  rot_and_project_safe <- function(i) {
    # Comprobar que pose y rotación son finitas
    if (!all(is.finite(c(ETdata$rX[i], ETdata$rY[i], ETdata$rZ[i],
                         ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])))) {
      return(c(NA_real_, NA_real_, NA_real_))
    }
    R <- Rcompose(ETdata$rX[i], ETdata$rY[i], ETdata$rZ[i])
    
    Li <- c(ETdata$goLx[i], ETdata$goLy[i], ETdata$goLz[i])
    Ri <- c(ETdata$goRx[i], ETdata$goRy[i], ETdata$goRz[i])
    Ld <- c(ETdata$gdLx[i], ETdata$gdLy[i], ETdata$gdLz[i])
    Rd <- c(ETdata$gdRx[i], ETdata$gdRy[i], ETdata$gdRz[i])
    
    if (!all(is.finite(c(Li, Ri, Ld, Rd)))) {
      return(c(NA_real_, NA_real_, NA_real_))
    }
    
    L_origin_w <- as.vector(R %*% Li) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
    R_origin_w <- as.vector(R %*% Ri) + c(ETdata$pX[i], ETdata$pY[i], ETdata$pZ[i])
    L_dir_w    <- as.vector(R %*% Ld)
    R_dir_w    <- as.vector(R %*% Rd)
    
    # pesos válidos (ya calculados fuera: Lw, Rw)
    if (!is.finite(Lw[i] + Rw[i]) || (Lw[i] + Rw[i]) <= 0) {
      return(c(NA_real_, NA_real_, NA_real_))
    }
    C_origin <- Lw[i]*L_origin_w + Rw[i]*R_origin_w
    C_dir    <- Lw[i]*L_dir_w    + Rw[i]*R_dir_w
    
    nrm <- sqrt(sum(C_dir^2))
    if (!is.finite(nrm) || nrm < 1e-9) return(c(NA_real_, NA_real_, NA_real_))
    C_dir <- C_dir / nrm
    
    ray_plane_intersect_safe(C_origin, C_dir, n = c(0,0,-1), p0 = c(0,0,3.2))
  }
  
  
  proj <- vapply(seq_len(nrow(ETdata)), rot_and_project_safe, FUN.VALUE = numeric(3))
  proj <- t(proj)
  ETdata$planeIntersect_x <- proj[,1]
  ETdata$planeIntersect_y <- proj[,2]
  ETdata$planeIntersect_z <- proj[,3]
}

if (!all(c("planeIntersect_x","planeIntersect_y") %in% names(ETdata))) {
  stop("Could not create planeIntersect_x/y. Check raw column names.")
}

## =========================
## 5) Invalid samples + padding and exclusions
## =========================
ETdata <- ETdata %>%
  dplyr::mutate(
    invalidData = ifelse(planeIntersect_x < -SHELF_BOUND | planeIntersect_x >  SHELF_BOUND |
                           planeIntersect_y < -SHELF_BOUND | planeIntersect_y >  SHELF_BOUND, 1, 0),
    invalidData = tidyr::replace_na(invalidData, 1)
  ) %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::mutate(
    invalidWindow = slider::slide_int(invalidData, ~ max(.x), .before = PAD_BLINK_N, .after = PAD_BLINK_N)
  ) %>%
  dplyr::ungroup()

## Trial-level loss
trial_loss <- ETdata %>% dplyr::group_by(uniqueID) %>% dplyr::summarise(propInvalid = mean(invalidWindow==1))
exclude_trials <- trial_loss %>% dplyr::filter(propInvalid > TRIAL_LOSS_T) %>% dplyr::pull(uniqueID)
keep_beeps <- ETdata %>% dplyr::filter(uniqueID %in% exclude_trials, beepDetectionTrial==1) %>% dplyr::pull(uniqueID) %>% unique()
ETdata <- ETdata %>% dplyr::filter(!(uniqueID %in% setdiff(exclude_trials, keep_beeps)))

## Subject-level loss
trial_counts <- ETdata %>%
  dplyr::distinct(uniqueID, subjectNr, trialNr) %>%
  dplyr::count(subjectNr, name = "n_kept") %>%
  dplyr::mutate(propInvalidTrials = 1 - n_kept / MAX_TRIALS)
exclude_subs <- trial_counts %>% dplyr::filter(propInvalidTrials > SUBJ_LOSS_T) %>% dplyr::pull(subjectNr)
ETdata <- ETdata %>% dplyr::filter(!(subjectNr %in% exclude_subs))

## Low-score trials (except beep trials)
scores <- ETdata %>% dplyr::group_by(uniqueID) %>% dplyr::slice(1) %>% dplyr::select(uniqueID, score, beepDetectionTrial)
low_score_trials <- scores %>% dplyr::filter(score <= SCORE_MIN, beepDetectionTrial==0) %>% dplyr::pull(uniqueID)
ETdata <- ETdata %>% dplyr::filter(!(uniqueID %in% low_score_trials))

## =========================
## 6) Velocities (raw, non-uniform time)
## =========================
vel_raw <- ETdata %>%
  dplyr::filter(invalidWindow==0) %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::mutate(
    dt_ms  = c(0, diff(t)),
    xdiff  = c(0, diff(planeIntersect_x)),
    ydiff  = c(0, diff(planeIntersect_y)),
    eucDist= sqrt(xdiff^2 + ydiff^2),
    euc_vel= eucDist / (dt_ms/1000),           # meters/second
    theta2D= rad2deg(atan2(ydiff, xdiff)),
    ang_vel= theta2D / (dt_ms/1000)            # degrees/second
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(uniqueID, t, dt_ms, euc_vel, ang_vel)

ETdata_clean <- ETdata %>%
  dplyr::select(-dt) %>%
  dplyr::left_join(vel_raw, by = c("uniqueID","t")) %>%
  dplyr::mutate(t_targ1StartTime = t - targ1StartTime) %>%
  dplyr::select(order(colnames(.)))

## =========================
## 7) Interpolation to a regular 120 Hz grid
## =========================
interp_one <- function(df){
  xi <- seq(0, max(df$t, na.rm=TRUE), by = DT_MS)
  px <- pracma::interp1(df$t, df$planeIntersect_x, xi, method = "linear")
  py <- pracma::interp1(df$t, df$planeIntersect_y, xi, method = "linear")
  out <- tibble::tibble(t = xi, planeIntersect_x = px, planeIntersect_y = py)
  out$subjectNr <- df$subjectNr[1]
  out$trialNr   <- df$trialNr[1]
  out$uniqueID  <- df$uniqueID[1]
  out
}

Interpdata <- ETdata_clean %>%
  dplyr::distinct(uniqueID, subjectNr, trialNr) %>%
  dplyr::pull(uniqueID) %>%
  purrr::map_dfr(\(id){
    cur <- ETdata_clean %>% dplyr::filter(uniqueID==id)
    suppressWarnings(interp_one(cur))
  })

## Add anchors and fixed dt; compute velocities on uniform grid
meta_tbl <- ETdata_clean %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(
    targ1StartTime   = dplyr::first(targ1StartTime),
    trialFinishTime  = dplyr::first(trialFinishTime),
    beepStartTime    = dplyr::first(beepStartTime),
    beepDetectTime   = dplyr::first(beepDetectTime),
    trialID          = dplyr::first(trialID),
    target0          = dplyr::first(target0),
    target1          = dplyr::first(target1),
    totalMoves       = dplyr::first(totalMoves),
    correctMoves     = dplyr::first(correctMoves),
    completion       = dplyr::first(completion),
    score            = dplyr::first(score),
    beepDetectionTrial = dplyr::first(beepDetectionTrial),
    .groups = "drop"
  )

Interpdata <- Interpdata %>%
  dplyr::left_join(meta_tbl, by = "uniqueID") %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::mutate(
    t_targ1StartTime     = t - targ1StartTime,
    t_beepStartTime      = t - beepStartTime,
    t_beepDetectionTime  = t - beepDetectTime,
    dt = DT_MS
  ) %>%
  dplyr::mutate(
    xdiff  = c(0, diff(planeIntersect_x)),
    ydiff  = c(0, diff(planeIntersect_y)),
    eucDist= sqrt(xdiff^2 + ydiff^2),
    euc_vel= eucDist / (dt/1000),        # m/s
    theta2D= rad2deg(atan2(ydiff, xdiff)),
    ang_vel= theta2D / (dt/1000)         # deg/s
  ) %>%
  dplyr::select(-xdiff, -ydiff, -eucDist, -theta2D) %>%
  dplyr::ungroup() %>%
  dplyr::select(order(colnames(.)))

## =========================
## 8) Save outputs + a tiny exclusion log
## =========================
save(ETdata_clean, file = file.path(DATA_PROCESSED, "ETdata_clean.Rdata"))
save(Interpdata,  file = file.path(DATA_PROCESSED, "Interpdata_clean.Rdata"))
readr::write_csv(ETdata_clean, file = file.path(DATA_PROCESSED, "ETdata_clean.csv"))
readr::write_csv(Interpdata,  file = file.path(DATA_PROCESSED, "Interpdata_clean.csv"))

exclusionTally <- list(
  trials_dropped_nonbeep = length(setdiff(exclude_trials, keep_beeps)),
  trials_kept_beep       = length(keep_beeps),
  subjects_dropped       = length(exclude_subs),
  trials_dropped_lowscore= length(low_score_trials)
)
writeLines(capture.output(str(exclusionTally)), con = file.path(RESULTS_DIR, "exclusion_tally.txt"))

message("✅ Preprocessing finished. Outputs in data/processed/")
