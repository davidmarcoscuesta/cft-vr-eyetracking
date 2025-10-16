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

## =========================
## 1) Load Trialdata
## =========================
Trialdata <- purrr::map_dfr(SUBS, \(s){
  fread(file.path("Data","IndividualSubjectData", sprintf("Trialdata_%d.txt", s)))
})
Trialdata <- Trialdata %>%
  filter(practiceNr == 3) %>%
  mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>%
  select(-practiceNr)

## Parse timestamps to clock::naive_time
Trialdata <- Trialdata %>%
  mutate(
    ymd = date_parse(str_extract(trialStartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    trialStartTime  = as_naive_time(year_month_day(get_year(ymd),get_month(ymd),get_day(ymd),
                                                   trialStartTime_h, trialStartTime_m, trialStartTime_s, trialStartTime_ms,
                                                   subsecond_precision = "millisecond")),
    ymd = date_parse(str_extract(targ1StartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    targ1StartTime  = as_naive_time(year_month_day(get_year(ymd),get_month(ymd),get_day(ymd),
                                                   targ1StartTime_h, targ1StartTime_m, targ1StartTime_s, targ1StartTime_ms,
                                                   subsecond_precision = "millisecond")),
    ymd = date_parse(str_extract(trialFinishTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    trialFinishTime = as_naive_time(year_month_day(get_year(ymd),get_month(ymd),get_day(ymd),
                                                   trialFinishTime_h, trialFinishTime_m, trialFinishTime_s, trialFinishTime_ms,
                                                   subsecond_precision = "millisecond"))
  ) %>%
  select(-matches("_(h|m|s|ms)$"), -ymd) %>%
  select(order(colnames(.)))

## =========================
## 2) Load ETdata
## =========================
# Option: load already-imported raw ETdata to save time
load(file.path("Data","rawETdata.Rdata"))  # creates ETdata

## Recompute trial-relative time (t) and inter-sample dt (ms)
ETdata <- ETdata %>%
  mutate(
    ymd = date_parse(str_extract(SysTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
    SysTime = as_naive_time(year_month_day(get_year(ymd),get_month(ymd),get_day(ymd),
                                           SysTime_h, SysTime_m, SysTime_s, SysTime_ms,
                                           subsecond_precision = "millisecond"))
  ) %>%
  group_by(uniqueID) %>%
  mutate(t = as.numeric(SysTime - SysTime[1]),
         dt = c(0, diff(t))) %>%
  ungroup() %>%
  select(-matches("SysTime_(h|m|s|ms)$"), -ymd, -time_stamp) %>%
  select(order(colnames(.)))

## Remove the single largest dt glitch and re-compute (empirically observed)
rm_row <- which.max(ETdata$dt) - 1
if(length(rm_row) == 1 && rm_row > 0){
  ETdata <- ETdata[-rm_row,] %>%
    group_by(uniqueID) %>%
    mutate(t = as.numeric(SysTime - SysTime[1]), dt = c(0, diff(t))) %>%
    ungroup()
}
## Known tail glitch in one trial
ETdata <- ETdata %>% filter(!(uniqueID == "8_7" & t >= 33663))

## Merge Trialdata
ETdata <- left_join(ETdata, Trialdata) %>% select(order(colnames(.)))

## Make trial-relative anchors (ms from trial start)
ETdata <- ETdata %>%
  group_by(uniqueID) %>%
  mutate(
    targ1StartTime  = as.numeric(targ1StartTime  - SysTime[1]),
    trialFinishTime = as.numeric(trialFinishTime - SysTime[1])
  ) %>%
  ungroup() %>%
  select(order(colnames(.)))

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
# If preProcessedETdata.Rdata (with planeIntersect_x/y) exists, load it:
if (file.exists(file.path("Data","preProcessedETdata.Rdata"))) {
  load(file.path("Data","preProcessedETdata.Rdata"))  # overwrites ETdata with planeIntersect_* columns
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
## 8) Save outputs + a tiny exclusion log
## =========================
dir.create("Data", showWarnings = FALSE, recursive = TRUE)
dir.create("Results", showWarnings = FALSE, recursive = TRUE)

save(ETdata_clean, file = file.path("Data","ETdata_clean.Rdata"))
save(Interpdata,  file = file.path("Data","Interpdata_clean.Rdata"))
write.csv(ETdata_clean, file = file.path("Data","ETdata_clean.csv"), row.names = FALSE)
write.csv(Interpdata,  file = file.path("Data","Interpdata_clean.csv"), row.names = FALSE)

exclusionTally <- list(
  trials_dropped_nonbeep = length(setdiff(exclude_trials, keep_beeps)),
  trials_kept_beep       = length(keep_beeps),
  subjects_dropped       = length(exclude_subs),
  trials_dropped_lowscore= length(low_score_trials)
)
writeLines(capture.output(str(exclusionTally)), con = file.path("Results","exclusion_tally.txt"))

message("✅ Preprocessing finished. Outputs in Data/ETdata_clean* and Data/Interpdata_clean*")
