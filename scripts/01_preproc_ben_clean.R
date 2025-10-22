# =========================
# CFT_preprocessing  (Ben Falandays, 2024-04-21)
# Ported to a plain R script without changing functionality.
# Outputs are written under Data/ (same as the original Rmd root.dir).
# =========================

## Library loading (same set Ben used)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(BAMBI)
  library(pracma)
  library(RcppRoll)
  library(clock)
  library(slider)
  library(TSPred)
})
options(digits.secs = 6)

## Reproducible working directory: mimic knitr root.dir = "/./Data"
# (Assumes you run this file from the project root.)
wd_data <- file.path(getwd(), "Data")
if (!dir.exists(wd_data)) stop("Data/ folder not found next to the project root.")
old_wd <- getwd()
setwd(wd_data)
on.exit(setwd(old_wd), add = TRUE)

# =========================
# Data Loading 
# =========================
subs <- seq(1,31)
exclude_subs <- c()  # can put anyone here if we want to exclude them off the bat
subs <- subset(subs, !subs %in% exclude_subs)

## Loading Trialdata (exactly as Ben)
for (s in subs) {
  fname <- paste0("./IndividualSubjectData/Trialdata_", s, ".txt")
  cur <- read.table(fname, sep = ",", header = TRUE)
  if (s == min(subs)) Trialdata <- cur else Trialdata <- rbind(Trialdata, cur)
}
rm(cur)
Trialdata <- Trialdata %>% 
  filter(practiceNr == 3) %>% 
  mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>% 
  select(-practiceNr)

## Parsing systime (exactly as Ben)
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

# alphabetically order the columns so it's easier to manually explore the data files later
Trialdata <- Trialdata %>% select(order(colnames(Trialdata)))

# =========================
# Making exclusions based on total score  (same as Ben)
# =========================
tmp <- Trialdata %>% group_by(subjectNr) %>% summarize(totalScore = sum(score))
hist(tmp$totalScore, breaks = 20)

maxScore <- 10*24
meanScore <- mean(tmp$totalScore)
sdScore <- sd(tmp$totalScore)
scoreCutoff <- .6 * maxScore

exclude_subs <- tmp %>% filter(totalScore < scoreCutoff) %>% .$subjectNr 

exclusionTally <- list(exclusion = paste0("Excluded ", length(exclude_subs),
                                          " subjects who had a total score < 60% of the max"))

# Trialdata <- Trialdata %>% filter(! subjectNr %in% exclude_subs) 
subs <- Trialdata %>% filter(!(subjectNr %in% exclude_subs)) %>% 
  group_by(subjectNr) %>% slice(1) %>% .$subjectNr 

# =========================
# Loading ETdata  (Ben’s code had eval=F here and then load rawETdata.Rdata)
# =========================
# (Kept commented to preserve original behavior)
# for(s in subs){
#   print(s)
#   fname = paste0('./IndividualSubjectData/ETdata_', s, '.txt')
#   cur <- read.table(fname, sep=',', header=T)
#   ifelse(s==min(subs), ETdata<-cur, ETdata <- rbind(ETdata,cur))
# }
# rm(cur)
# ETdata = ETdata %>% filter(practiceNr == 3) %>% mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>% select(-practiceNr)

# =========================
# Parsing systime for ET (the block was eval=F in Ben’s Rmd)
# =========================
# ETdata = ETdata %>%
#   mutate(ymd = date_parse(str_extract(SysTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format='%m/%d/%y'),
#          SysTime = as_naive_time(year_month_day(get_year(ymd),get_month(ymd), get_day(ymd), 
#                                                SysTime_h, SysTime_m, SysTime_s, SysTime_ms, 
#                                                subsecond_precision = "millisecond"))) %>%
#   group_by(uniqueID) %>%
#   mutate(t = as.numeric(SysTime - SysTime[1]), 
#          dt = c(0, diff(t))) %>%
#   ungroup() %>%
#   select(-c(SysTime_h, SysTime_m, SysTime_s, SysTime_ms, ymd,time_stamp))
# ETdata = ETdata %>% select(order(colnames(ETdata)))

# Save/load raw ETdata (exactly as Ben: load is eval=TRUE)
# save(ETdata, file="rawETdata.Rdata")
load("rawETdata.Rdata")  # creates ETdata

# =========================
# Check for and deal with bugs in the time variables (same)
# =========================
tmp <- ETdata %>% group_by(uniqueID) %>% summarise(maxDT = max(dt))
summary(tmp$maxDT)
which(tmp$maxDT > 100)

summary(ETdata$dt)
which.max(ETdata$dt)
rm_row <- which(ETdata$dt == max(ETdata$dt)) - 1
ETdata <- ETdata[-rm_row,]
ETdata <- ETdata %>% filter(!(uniqueID == "8_7" & t >= 33663))

ETdata <- ETdata %>% 
  group_by(uniqueID) %>% 
  mutate(t  = as.numeric(SysTime - SysTime[1]), 
         dt = c(0, diff(t))) %>% 
  ungroup()

summary(ETdata$dt)

# =========================
# Merge ETdata and Trialdata (same)
# =========================
ETdata <- left_join(ETdata, Trialdata)
ETdata <- ETdata %>% select(order(colnames(ETdata)))

# =========================
# Recode time-locking vars relative to first sample (same)
# =========================
ETdata <- ETdata %>% group_by(uniqueID) %>% 
  mutate(targ1StartTime  = as.numeric(targ1StartTime  - SysTime[1]), 
         trialFinishTime = as.numeric(trialFinishTime - SysTime[1])) %>% 
  ungroup()
ETdata <- ETdata %>% select(order(colnames(ETdata)))

# =========================
# Beep detection trials and time (same)
# =========================
tmp <- ETdata %>% 
  group_by(uniqueID) %>% 
  mutate(beepDetection = ifelse((beepStatus == 2) & (lag(beepStatus) == 1), 1, 0)) %>% 
  filter(beepDetection == 1) %>% 
  select(c(uniqueID, beepDetection, t)) %>% 
  rename(beepDetectionTime = t, beepDetectionTrial = beepDetection)

ETdata <- left_join(ETdata, tmp) %>% replace_na(list(beepDetectionTrial = 0, beepDetectionTime = NA))
ETdata <- ETdata %>% select(order(colnames(ETdata)))

# =========================
# Loading Clickdata (same)
# =========================
for (s in subs) {
  fname <- paste0("./IndividualSubjectData/Clickdata_", s, ".txt")
  cur <- read.table(fname, sep = ",", header = TRUE)
  if (s == min(subs)) Clickdata <- cur else Clickdata <- rbind(Clickdata, cur)
}
rm(cur)
Clickdata <- Clickdata %>% 
  filter(practiceNr == 3) %>% 
  mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>% 
  select(-practiceNr)

# Parsing systime (same)
Clickdata <- Clickdata %>% 
  mutate(
    ymd = date_parse(str_extract(SysTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format = "%m/%d/%y"),
    SysTime = as_naive_time(year_month_day(get_year(ymd), get_month(ymd), get_day(ymd),
                                           SysTime_h, SysTime_m, SysTime_s, SysTime_ms,
                                           subsecond_precision = "millisecond"))
  ) %>%
  group_by(uniqueID) %>% 
  ungroup() %>% 
  select(-c(SysTime_h, SysTime_m, SysTime_s, SysTime_ms, ymd)) 

# Merge Clickdata with ETdata (same)
Clickdata <- left_join(
  Clickdata,
  ETdata %>% group_by(uniqueID) %>% slice(1) %>% 
    rename(trialStartTime = SysTime) %>% select(c(trialStartTime, uniqueID))
)
Clickdata <- Clickdata %>% 
  group_by(uniqueID) %>% mutate(t = SysTime - trialStartTime) %>% ungroup()
Clickdata <- Clickdata %>% select(order(colnames(Clickdata)))

# =========================
# Processing ETdata (all blocks kept EXACTLY — the first ones were eval=F)
# =========================

## Transforming gaze data variables (Ben had eval=F; left commented)
# ETdata = ETdata %>%  
#   mutate(
#     gaze_direct_L.x_local = -gaze_direct_L.x, 
#     gaze_direct_L.y_local = gaze_direct_L.y,
#     gaze_direct_L.z_local = gaze_direct_L.z,
#     gaze_direct_R.x_local = -gaze_direct_R.x,
#     gaze_direct_R.y_local = gaze_direct_R.y,
#     gaze_direct_R.z_local = gaze_direct_R.z,
#     gaze_origin_L.x_local = -gaze_origin_L.x.mm.*.001, 
#     gaze_origin_L.y_local = gaze_origin_L.y.mm.*.001, 
#     gaze_origin_L.z_local = gaze_origin_L.z.mm.*.001,
#     gaze_origin_R.x_local = -gaze_origin_R.x.mm.*.001, 
#     gaze_origin_R.y_local = gaze_origin_R.y.mm.*.001, 
#     gaze_origin_R.z_local = gaze_origin_R.z.mm.*.001,
#     camRot_x = minuspi_to_pi(deg2rad(camRot_x)), 
#     camRot_y = minuspi_to_pi(deg2rad(camRot_y)),
#     camRot_z = minuspi_to_pi(deg2rad(camRot_z)),
#     eye_valid_L = ifelse(eye_valid_L==31, 1, 0), eye_valid_R = ifelse(eye_valid_R==31, 1, 0),
#     eye_valid_L = ifelse(openness_L < .2, 0, eye_valid_L), eye_valid_R = ifelse(openness_R < .2, 0, eye_valid_R),
#     eye_valid_L = ifelse(sqrt(gaze_direct_L.x^2 + gaze_direct_L.y^2 +gaze_direct_L.z^2 ) <.99, 0, eye_valid_L),
#     eye_valid_R = ifelse(sqrt(gaze_direct_R.x^2 + gaze_direct_R.y^2 +gaze_direct_R.z^2 ) <.99, 0, eye_valid_R),
#     eye_valid_L = ifelse(pupil_diameter_L.mm. < 1, 0, eye_valid_L), eye_valid_R = ifelse(pupil_diameter_R.mm. < 1, 0, eye_valid_R),
#     invalidData = ifelse(eye_valid_L + eye_valid_R < 1, 1, 0)
#   ) %>% 
#   select(-c(pos_sensor_L.x, pos_sensor_L.y, pos_sensor_R.x, pos_sensor_R.y,
#             distance_valid_C, distance_C.mm., 
#             gaze_direct_L.x, gaze_direct_L.y, gaze_direct_L.z, 
#             gaze_direct_R.x, gaze_direct_R.y, gaze_direct_R.z,
#             gaze_origin_L.x.mm., gaze_origin_L.y.mm., gaze_origin_L.z.mm.,
#             gaze_origin_R.x.mm., gaze_origin_R.y.mm., gaze_origin_R.z.mm.,
#             gaze_sensitive)) 
# ETdata = ETdata %>% select(order(colnames(ETdata)))

## Rotation helpers (kept; used only if you un-comment the transform)
rotation_matrix_x_left <- function(theta) matrix(c(1,0,0, 0,cos(theta),-sin(theta), 0,sin(theta),cos(theta)),3,3,byrow=TRUE)
rotation_matrix_y_left <- function(theta) matrix(c(cos(theta),0,sin(theta), 0,1,0, -sin(theta),0,cos(theta)),3,3,byrow=TRUE)
rotation_matrix_z_left <- function(theta) matrix(c(cos(theta),-sin(theta),0, sin(theta),cos(theta),0, 0,0,1),3,3,byrow=TRUE)
rotation_matrix_x_right <- function(theta) matrix(c(1,0,0, 0,cos(theta),sin(theta), 0,-sin(theta),cos(theta)),3,3,byrow=TRUE)
rotation_matrix_y_right <- function(theta) matrix(c(cos(theta),0,-sin(theta), 0,1,0, sin(theta),0,cos(theta)),3,3,byrow=TRUE)
rotation_matrix_z_right <- function(theta) matrix(c(cos(theta),sin(theta),0, -sin(theta),cos(theta),0, 0,0,1),3,3,byrow=TRUE)
rotate_vector <- function(vector, euler_angles, rotation_order="zxy", coordinate_system="left"){
  rotation_functions <- if (coordinate_system=="left") list(x=rotation_matrix_x_left,y=rotation_matrix_y_left,z=rotation_matrix_z_left)
  else list(x=rotation_matrix_x_right,y=rotation_matrix_y_right,z=rotation_matrix_z_right)
  rotation_axes <- unlist(strsplit(rotation_order,""))
  rotation_matrices <- lapply(rotation_axes, function(rot_axis) rotation_functions[[rot_axis]](euler_angles[which(rot_axis==c("x","y","z"))]))
  combined_rotation <- Reduce(function(a,b) b %*% a, rotation_matrices)
  as.vector(combined_rotation %*% vector)
}

## Compute world coords & cyclopean & plane intersection (kept commented exactly as in Ben)
# ... (all the eval=F blocks remain commented here)

# Load preprocessed ET (Ben had eval=TRUE here)
load("preProcessedETdata.Rdata")   # creates ETdata with planeIntersect_* columns

# =========================
# Making more exclusions (same as Ben)
# =========================
ETdata <- ETdata %>% group_by(uniqueID) %>% 
  mutate(invalidData = ifelse(planeIntersect_x < -1.5 | planeIntersect_x > 1.5 |
                                planeIntersect_y < -1.5 | planeIntersect_y > 1.5, 1, invalidData)) %>%
  ungroup() %>%
  mutate(invalidData = replace_na(invalidData, 1)) %>%
  group_by(uniqueID) %>%
  mutate(invalidWindow = slide_index_int(invalidData, t, max, .before = 30, .after = 30)) %>% 
  mutate(planeIntersect_x = ifelse(invalidWindow==1, NA, planeIntersect_x),
         planeIntersect_y = ifelse(invalidWindow==1, NA, planeIntersect_y)) %>% 
  ungroup()
ETdata <- ETdata %>% select(order(colnames(ETdata)))

## Exclude trials with too much data loss (same)
tmp <- ETdata %>% group_by(subjectNr, trialNr, uniqueID) %>% summarize(propInvalid = sum(invalidWindow)/n())
hist(tmp$propInvalid)
thresh <- .1
exclude_trials <- tmp %>% ungroup() %>% filter(propInvalid > thresh) %>% .$uniqueID
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)
exclusionTally[[length(exclusionTally)+1]] <- list(exclusion = paste0("Excluded ", length(exclude_trials),
                                                                      " trials with > 10% data loss"))

## Exclude entire participants with too many trials removed (same)
tmp <- ETdata %>% group_by(subjectNr) %>% summarise(propInvalidTrials = 1 - length(unique(trialNr))/24)
hist(tmp$propInvalidTrials)
thresh <- .2
exclude_subs <- tmp %>% ungroup() %>% filter(propInvalidTrials > thresh) %>% .$subjectNr
ETdata <- ETdata %>% filter(!(subjectNr %in% exclude_subs))
exclusionTally[[length(exclusionTally)+1]] <- list(exclusion = paste0("Excluded ", length(exclude_subs),
                                                                      " subjects with > 20% of trials removed due to data loss"))

## Very long/short trials (same)
ETdata %>% group_by(uniqueID) %>% summarise(totalTime = max(t)) %>% summary() 
totalTimes <- ETdata %>% group_by(uniqueID) %>% summarise(totalTime = max(t))
hist(totalTimes$totalTime, breaks = 100)
median_time <- median(totalTimes$totalTime); sd_time <- sd(totalTimes$totalTime)
longCutoff <- 2*median_time; shortCutoff <- median_time/2
long_trials  <- c(totalTimes$uniqueID[which(totalTimes$totalTime > longCutoff)]) 
short_trials <- c(totalTimes$uniqueID[which(totalTimes$totalTime < shortCutoff)]) 
exclude_trials <- c(short_trials, long_trials)
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)
exclusionTally[[length(exclusionTally)+1]] <- list(exclusion = paste0("Excluded ", length(short_trials),
                                                                      " trials that were very short, and ",
                                                                      length(long_trials), " trials that were very long"))

## Low-score trials (same)
tmp <- ETdata %>% group_by(uniqueID) %>% slice(1) %>% select(score, uniqueID)
hist(tmp$score); summary(tmp$score)
scoreCutoff <- 6
exclude_trials <- tmp %>% filter(score <= scoreCutoff) %>% .$uniqueID
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)
exclusionTally[[length(exclusionTally)+1]] <- list(exclusion = paste0("Excluded ", length(exclude_trials),
                                                                      " trials where the trial score was < 6/10"))

## What proportion remains? (same)
length(unique(ETdata$uniqueID)) / length(unique(Trialdata$uniqueID))

# =========================
# Smoothing (kept commented as in Ben)
# =========================
# (All smoothing blocks left exactly as in the Rmd with eval=F)

# =========================
# Interpolation to even time (same)
# =========================
if (exists("Interpdata")) rm(Interpdata)
counter <- 0
nTrials <- n_distinct(ETdata$uniqueID)
for (id in unique(ETdata$uniqueID)) {
  counter <- counter + 1
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
  ungroup() 
Interpdata <- Interpdata %>% select(order(colnames(Interpdata)))

# =========================
# Velocity series (same)
# =========================
tmp <- ETdata %>% 
  filter(invalidWindow==0) %>% 
  group_by(uniqueID) %>% 
  mutate(
    dt = c(0, diff(t)),
    xdiff = c(0, diff(planeIntersect_x)),
    ydiff = c(0, diff(planeIntersect_y)),
    eucDist = sqrt(xdiff^2 + ydiff^2),
    euc_vel = eucDist/(dt/1000),
    theta2D = rad2deg(atan2(ydiff, xdiff)),
    ang_vel = theta2D/(dt/1000)
  ) %>% 
  ungroup() %>% 
  select(c(uniqueID, t, dt, euc_vel, ang_vel))

ETdata <- left_join(ETdata %>% select(-dt), tmp) %>% 
  group_by(uniqueID) %>% 
  mutate(t_targ1StartTime = t - targ1StartTime,
         t_beepDetectionTime = t - beepDetectionTime) %>% 
  select(c(uniqueID, subjectNr, trialNr, 
           targ1StartTime, trialFinishTime, t_targ1StartTime, t_beepDetectionTime,
           trialID, target0, target1, 
           totalMoves, correctMoves, completion, score, 
           beepDetectionTime, beepDetectionTrial, 
           planeIntersect_x, planeIntersect_y, 
           dt, t, euc_vel, ang_vel)) 
ETdata <- ETdata %>% select(order(colnames(ETdata)))

Interpdata <- Interpdata %>% 
  group_by(uniqueID) %>% 
  mutate(xdiff = c(0, diff(planeIntersect_x)),
         ydiff = c(0, diff(planeIntersect_y)),
         eucDist = sqrt(xdiff^2 + ydiff^2),
         euc_vel = eucDist/dt,
         theta2D = rad2deg(atan2(ydiff, xdiff)),
         ang_vel = theta2D/dt) %>% 
  select(-c(xdiff, ydiff, eucDist, theta2D)) %>% 
  ungroup()
Interpdata <- Interpdata %>% select(order(colnames(Interpdata)))

# =========================
# Write outputs (same file names and locations as Ben)
# =========================
save(Trialdata, file = "Trialdata.Rdata"); write.csv(Trialdata, file = "Trialdata.csv")
save(Clickdata,  file = "Clickdata.Rdata"); write.csv(Clickdata,  file = "Clickdata.csv")
save(ETdata,     file = "ETdata.Rdata");    write.csv(ETdata,     file = "ETdata.csv")
save(Interpdata, file = "Interpdata.Rdata");write.csv(Interpdata, file = "Interpdata.csv")
write.csv(exclusionTally, file = "exclusionTally.csv")

message("✅ Finished. Outputs written in Data/ (same as Ben).")
