# =========================
# CFT_preprocessing  (Ben Falandays, 2024-04-21)
# Plain R script adapted to our folders: data/raw, data/processed, results.
# Functionality, variable names and logic are unchanged.
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(BAMBI)
  library(pracma)
  library(RcppRoll)
  library(clock)
  options(digits.secs = 6)
  library(slider)
  library(TSPred)
})
options(digits.secs = 6)

# ---- Paths (only change vs Ben) ----
DATA_DIR       <- "data"
RAW_DIR        <- file.path(DATA_DIR, "raw")
RAW_ALT_DIR    <- file.path(RAW_DIR, "IndividualSubjectData")  # fallback if files are there
PROC_DIR       <- file.path(DATA_DIR, "processed")
RESULTS_DIR    <- "results"
dir.create(PROC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Helper to find a file in raw root or in IndividualSubjectData/
raw_path <- function(filename) {
  p1 <- file.path(RAW_DIR, filename)
  if (file.exists(p1)) return(p1)
  p2 <- file.path(RAW_ALT_DIR, filename)
  if (file.exists(p2)) return(p2)
  return(p1) # default
}

# =========================
# Data Loading 
# =========================
subs <- seq(1,31)
exclude_subs <- c()
subs <- subset(subs, !subs %in% exclude_subs)

## Loading Trialdata (same logic as Ben; just paths)
for (s in subs) {
  fname <- raw_path(paste0("Trialdata_", s, ".txt"))
  cur <- read.table(fname, sep = ",", header = TRUE)
  if (s == min(subs)) Trialdata <- cur else Trialdata <- rbind(Trialdata, cur)
}
rm(cur)
Trialdata <- Trialdata %>% 
  filter(practiceNr == 3) %>% 
  mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>% 
  select(-practiceNr)

## Parsing systime (unchanged)
Trialdata <- Trialdata %>% 
  mutate(
    ymd = date_parse(str_extract(trialStartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    trialStartTime = as_naive_time(year_month_day(get_year(ymd), get_month(ymd), get_day(ymd),
                                                  trialStartTime_h, trialStartTime_m, trialStartTime_s, trialStartTime_ms,
                                                  subsecond_precision = "millisecond"))
  ) %>% 
  mutate(
    ymd = date_parse(str_extract(targ1StartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    targ1StartTime = as_naive_time(year_month_day(get_year(ymd), get_month(ymd), get_day(ymd),
                                                  targ1StartTime_h, targ1StartTime_m, targ1StartTime_s, targ1StartTime_ms,
                                                  subsecond_precision = "millisecond"))
  ) %>% 
  mutate(
    ymd = date_parse(str_extract(trialFinishTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    trialFinishTime = as_naive_time(year_month_day(get_year(ymd), get_month(ymd), get_day(ymd),
                                                   trialFinishTime_h, trialFinishTime_m, trialFinishTime_s, trialFinishTime_ms,
                                                   subsecond_precision = "millisecond"))
  ) %>% 
  select(-c(trialStartTime_h, trialStartTime_m, trialStartTime_s, trialStartTime_ms, trialStartTime,
            targ1StartTime_h, targ1StartTime_m, targ1StartTime_s, targ1StartTime_ms, 
            trialFinishTime_h, trialFinishTime_m, trialFinishTime_s, trialFinishTime_ms, ymd))

Trialdata <- Trialdata %>% select(order(colnames(Trialdata)))

# =========================
# Exclusions based on total score (unchanged)
# =========================
tmp <- Trialdata %>% group_by(subjectNr) %>% summarize(totalScore = sum(score))
hist(tmp$totalScore, breaks = 20)

maxScore <- 10*24
scoreCutoff <- .6 * maxScore
exclude_subs <- tmp %>% filter(totalScore < scoreCutoff) %>% .$subjectNr 

exclusionTally <- list(exclusion = paste0("Excluded ", length(exclude_subs),
                                          " subjects who had a total score < 60% of the max"))

subs <- Trialdata %>% filter(!(subjectNr %in% exclude_subs)) %>% 
  group_by(subjectNr) %>% slice(1) %>% .$subjectNr 

# =========================
# ETdata (keep Ben’s behavior: load pre-saved)
# =========================
# The Rmd had the import loop commented out and then: load("rawETdata.Rdata")
load(raw_path("rawETdata.Rdata"))  # creates ETdata

# Check and fix time glitches (unchanged)
tmp <- ETdata %>% group_by(uniqueID) %>% summarise(maxDT = max(dt))
summary(tmp$maxDT); which(tmp$maxDT > 100)
summary(ETdata$dt); which.max(ETdata$dt)
rm_row <- which(ETdata$dt == max(ETdata$dt)) - 1
ETdata <- ETdata[-rm_row,]
ETdata <- ETdata %>% filter(!(uniqueID == "8_7" & t >= 33663))
ETdata <- ETdata %>% 
  group_by(uniqueID) %>% 
  mutate(t = as.numeric(SysTime - SysTime[1]), dt = c(0, diff(t))) %>% 
  ungroup()
summary(ETdata$dt)

# Merge with Trialdata; recode anchors (unchanged)
ETdata <- left_join(ETdata, Trialdata) %>% select(order(colnames(.)))
ETdata <- ETdata %>% group_by(uniqueID) %>% 
  mutate(targ1StartTime = as.numeric(targ1StartTime - SysTime[1]), 
         trialFinishTime = as.numeric(trialFinishTime - SysTime[1])) %>% 
  ungroup() %>% select(order(colnames(.)))

# Beep detection (unchanged)
tmp <- ETdata %>% 
  group_by(uniqueID) %>% 
  mutate(beepDetection = ifelse((beepStatus == 2) & (lag(beepStatus) == 1), 1, 0)) %>% 
  filter(beepDetection == 1) %>% 
  select(c(uniqueID, beepDetection, t)) %>% 
  rename(beepDetectionTime = t, beepDetectionTrial = beepDetection)
ETdata <- left_join(ETdata, tmp) %>% replace_na(list(beepDetectionTrial = 0, beepDetectionTime = NA))
ETdata <- ETdata %>% select(order(colnames(.)))

# =========================
# Clickdata (unchanged, just paths)
# =========================
for (s in subs) {
  fname <- raw_path(paste0("Clickdata_", s, ".txt"))
  cur <- read.table(fname, sep = ",", header = TRUE)
  if (s == min(subs)) Clickdata <- cur else Clickdata <- rbind(Clickdata, cur)
}
rm(cur)
Clickdata <- Clickdata %>% 
  filter(practiceNr == 3) %>% 
  mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>% select(-practiceNr)

Clickdata <- Clickdata %>% 
  mutate(
    ymd = date_parse(str_extract(SysTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    SysTime = as_naive_time(year_month_day(get_year(ymd), get_month(ymd), get_day(ymd),
                                           SysTime_h, SysTime_m, SysTime_s, SysTime_ms,
                                           subsecond_precision = "millisecond"))
  ) %>% select(-c(SysTime_h, SysTime_m, SysTime_s, SysTime_ms, ymd)) 

Clickdata <- left_join(
  Clickdata,
  ETdata %>% group_by(uniqueID) %>% slice(1) %>% 
    rename(trialStartTime = SysTime) %>% select(c(trialStartTime, uniqueID))
) %>% 
  group_by(uniqueID) %>% mutate(t = SysTime - trialStartTime) %>% ungroup() %>% 
  select(order(colnames(.)))

# =========================
# Plane projection & preprocessing (Ben loads preProcessedETdata.Rdata)
# =========================
load(raw_path("preProcessedETdata.Rdata"))  # creates ETdata with planeIntersect_*

# =========================
# Exclusions (unchanged)
# =========================
ETdata <- ETdata %>% group_by(uniqueID) %>% 
  mutate(invalidData = ifelse(planeIntersect_x < -1.5 | planeIntersect_x >  1.5 |
                                planeIntersect_y < -1.5 | planeIntersect_y >  1.5, 1, invalidData)) %>%
  ungroup() %>% mutate(invalidData = replace_na(invalidData, 1)) %>%
  group_by(uniqueID) %>% 
  mutate(invalidWindow = slide_index_int(invalidData, t, max, .before = 30, .after = 30),
         planeIntersect_x = ifelse(invalidWindow==1, NA, planeIntersect_x),
         planeIntersect_y = ifelse(invalidWindow==1, NA, planeIntersect_y)) %>%
  ungroup() %>% select(order(colnames(.)))

tmp <- ETdata %>% group_by(subjectNr, trialNr, uniqueID) %>% 
  summarize(propInvalid = sum(invalidWindow)/n())
hist(tmp$propInvalid)
thresh <- .1
exclude_trials <- tmp %>% ungroup() %>% filter(propInvalid > thresh) %>% .$uniqueID
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)

exclusionTally[[length(exclusionTally)+1]] <- list(
  exclusion = paste0("Excluded ", length(exclude_trials), " trials with > 10% data loss")
)

tmp <- ETdata %>% group_by(subjectNr) %>% summarise(propInvalidTrials = 1 - length(unique(trialNr))/24)
hist(tmp$propInvalidTrials)
thresh <- .2
exclude_subs <- tmp %>% ungroup() %>% filter(propInvalidTrials > thresh) %>% .$subjectNr
ETdata <- ETdata %>% filter(!(subjectNr %in% exclude_subs))
exclusionTally[[length(exclusionTally)+1]] <- list(
  exclusion = paste0("Excluded ", length(exclude_subs), " subjects with > 20% of trials removed due to data loss")
)

ETdata %>% group_by(uniqueID) %>% summarise(totalTime = max(t)) %>% summary()
totalTimes <- ETdata %>% group_by(uniqueID) %>% summarise(totalTime = max(t))
hist(totalTimes$totalTime, breaks = 100)
median_time <- median(totalTimes$totalTime)
longCutoff <- 2*median_time; shortCutoff <- median_time/2
long_trials  <- totalTimes$uniqueID[which(totalTimes$totalTime > longCutoff)] 
short_trials <- totalTimes$uniqueID[which(totalTimes$totalTime < shortCutoff)] 
exclude_trials <- c(short_trials, long_trials)
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)
exclusionTally[[length(exclusionTally)+1]] <- list(
  exclusion = paste0("Excluded ", length(short_trials), " trials that were very short, and ",
                     length(long_trials), " trials that were very long")
)

tmp <- ETdata %>% group_by(uniqueID) %>% slice(1) %>% select(score, uniqueID)
hist(tmp$score); summary(tmp$score)
scoreCutoff <- 6
exclude_trials <- tmp %>% filter(score <= scoreCutoff) %>% .$uniqueID
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)
exclusionTally[[length(exclusionTally)+1]] <- list(
  exclusion = paste0("Excluded ", length(exclude_trials), " trials where the trial score was < 6/10")
)

length(unique(ETdata$uniqueID)) / length(unique(Trialdata$uniqueID))

# =========================
# Interpolation (unchanged)
# =========================
if (exists("Interpdata")) rm(Interpdata)
for (id in unique(ETdata$uniqueID)) {
  cur <- ETdata %>% filter(uniqueID == id)
  interp_points <- seq(0, max(cur$t), 1000/120)
  planeIntersect_x <- pracma::interp1(x = cur$t, y = cur$planeIntersect_x, xi = interp_points, method = "linear")
  planeIntersect_y <- pracma::interp1(x = cur$t, y = cur$planeIntersect_y, xi = interp_points, method = "linear")
  newdata <- data.frame(t = interp_points, planeIntersect_x = planeIntersect_x, planeIntersect_y = planeIntersect_y)
  newdata$subjectNr <- unique(cur$subjectNr); newdata$trialNr <- unique(cur$trialNr); newdata$uniqueID <- unique(cur$uniqueID)
  if (!exists("Interpdata")) Interpdata <- newdata else Interpdata <- rbind(Interpdata, newdata)
}
Interpdata <- left_join(
  Interpdata, 
  ETdata %>% group_by(uniqueID) %>% 
    select(c(uniqueID, targ1StartTime, trialFinishTime, trialID, target0, target1,
             totalMoves, correctMoves, completion, score, beepDetectionTime, beepDetectionTrial)) %>% 
    slice(1)
)
Interpdata <- Interpdata %>% group_by(uniqueID) %>% 
  mutate(t_targ1StartTime = t - targ1StartTime,
         t_beepDetectionTime = t - beepDetectionTime,
         dt = c(0, diff(t))) %>% 
  ungroup() %>% 
  select(order(colnames(.)))

Interpdata <- Interpdata %>% 
  group_by(uniqueID) %>% 
  mutate(xdiff = c(0, diff(planeIntersect_x)),
         ydiff = c(0, diff(planeIntersect_y)),
         eucDist = sqrt(xdiff^2 + ydiff^2),
         euc_vel = eucDist/dt,
         theta2D = rad2deg(atan2(ydiff, xdiff)),
         ang_vel = theta2D/dt)) %>% 
  select(-c(xdiff, ydiff, eucDist, theta2D)) %>% 
  ungroup() %>% 
  select(order(colnames(.)))

# =========================
# Save outputs (same names; now under data/processed)
# =========================
save(Trialdata, file = file.path(PROC_DIR, "Trialdata.Rdata"))
write.csv(Trialdata, file.path(PROC_DIR, "Trialdata.csv"), row.names = FALSE)

save(Clickdata,  file = file.path(PROC_DIR, "Clickdata.Rdata"))
write.csv(Clickdata, file.path(PROC_DIR, "Clickdata.csv"), row.names = FALSE)

save(ETdata,     file = file.path(PROC_DIR, "ETdata.Rdata"))
write.csv(ETdata, file.path(PROC_DIR, "ETdata.csv"), row.names = FALSE)

save(Interpdata, file = file.path(PROC_DIR, "Interpdata.Rdata"))
write.csv(Interpdata, file.path(PROC_DIR, "Interpdata.csv"), row.names = FALSE)

write.csv(exclusionTally, file = file.path(PROC_DIR, "exclusionTally.csv"), row.names = FALSE)

message("✅ Finished. Inputs read from data/raw, outputs written to data/processed.")
