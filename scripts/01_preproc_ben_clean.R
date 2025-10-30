# =========================
# CFT_preprocessing (Ben) – rutas adaptadas a data/raw & data/processed
# Funcionalidad intacta (mismos nombres y lógica que el Rmd de Ben)
# =========================

suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(dplyr)
  library(ggplot2); library(BAMBI); library(pracma)
  library(RcppRoll); library(clock); library(slider); library(TSPred)
})
options(digits.secs = 6)

# ---- Directorios ----
RAW_DIR     <- "data/raw"
PROC_DIR    <- "data/processed"
dir.create(PROC_DIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# Data Loading 
# =========================
subs <- seq(1,31)
exclude_subs <- c()
subs <- subset(subs, !subs %in% exclude_subs)

## Trialdata (idéntico a Ben, cambiando solo la ruta)
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

## Parsing systime (igual que Ben)
Trialdata <- Trialdata %>% 
  mutate(ymd = date_parse(str_extract(trialStartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
         trialStartTime = as_naive_time(year_month_day(get_year(ymd),get_month(ymd),get_day(ymd),
                                                       trialStartTime_h, trialStartTime_m, trialStartTime_s, trialStartTime_ms,
                                                       subsecond_precision = "millisecond"))) %>% 
  mutate(ymd = date_parse(str_extract(targ1StartTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
         targ1StartTime = as_naive_time(year_month_day(get_year(ymd),get_month(ymd),get_day(ymd),
                                                       targ1StartTime_h, targ1StartTime_m, targ1StartTime_s, targ1StartTime_ms,
                                                       subsecond_precision = "millisecond"))) %>% 
  mutate(ymd = date_parse(str_extract(trialFinishTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
         trialFinishTime = as_naive_time(year_month_day(get_year(ymd),get_month(ymd),get_day(ymd),
                                                        trialFinishTime_h, trialFinishTime_m, trialFinishTime_s, trialFinishTime_ms,
                                                        subsecond_precision = "millisecond"))) %>% 
  select(-c(trialStartTime_h, trialStartTime_m,trialStartTime_s, trialStartTime_ms, trialStartTime,
            targ1StartTime_h, targ1StartTime_m,targ1StartTime_s, targ1StartTime_ms, 
            trialFinishTime_h, trialFinishTime_m,trialFinishTime_s, trialFinishTime_ms, ymd)) %>% 
  select(order(colnames(.)))

# =========================
# Exclusión por score total (igual que Ben)
# =========================
tmp <- Trialdata %>% group_by(subjectNr) %>% summarize(totalScore = sum(score))
maxScore <- 10*24
scoreCutoff <- .6*maxScore
exclude_subs <- tmp %>% filter(totalScore < scoreCutoff) %>% .$subjectNr
exclusionTally <- list(exclusion=paste0("Excluded ", length(exclude_subs),
                                        " subjects who had a total score < 60% of the max"))
subs <- Trialdata %>% filter(!(subjectNr %in% exclude_subs)) %>% group_by(subjectNr) %>% slice(1) %>% .$subjectNr 

# =========================
# ETdata
# =========================
# El Rmd de Ben tenía el bucle de importación comentado y luego:
#   save(ETdata,"rawETdata.Rdata")  (opcional)  y después:
#   load("rawETdata.Rdata")
load(raw_path("rawETdata.Rdata"))  # crea ETdata

# Glitches de tiempo (igual que Ben)
tmp <- ETdata %>% group_by(uniqueID) %>% summarise(maxDT = max(dt))
summary(ETdata$dt); which.max(ETdata$dt)
rm_row <- which(ETdata$dt==max(ETdata$dt)) - 1
ETdata <- ETdata[-rm_row,]
ETdata <- ETdata %>% filter(!(uniqueID == "8_7" & t >= 33663)) %>% 
  group_by(uniqueID) %>% 
  mutate(t = as.numeric(SysTime - SysTime[1]), dt = c(0, diff(t))) %>% 
  ungroup()

# Merge con Trialdata + anclajes relativos (igual)
ETdata <- left_join(ETdata, Trialdata) %>% select(order(colnames(.)))
ETdata <- ETdata %>% group_by(uniqueID) %>% 
  mutate(targ1StartTime = as.numeric(targ1StartTime - SysTime[1]), 
         trialFinishTime = as.numeric(trialFinishTime - SysTime[1])) %>% 
  ungroup() %>% select(order(colnames(.)))

# Beeps (igual)
tmp <- ETdata %>% 
  group_by(uniqueID) %>% 
  mutate(beepDetection = ifelse((beepStatus == 2 ) & (lag(beepStatus) == 1), 1, 0)) %>% 
  filter(beepDetection == 1) %>% 
  select(c(uniqueID, beepDetection, t)) %>% 
  rename(beepDetectionTime = t, beepDetectionTrial = beepDetection)
ETdata <- left_join(ETdata, tmp) %>% replace_na(list(beepDetectionTrial = 0, beepDetectionTime = NA)) %>% 
  select(order(colnames(.)))

# =========================
# Clickdata (igual; solo cambia path)
# =========================
for (s in subs) {
  fname <- raw_path(paste0("Clickdata_", s, ".txt"))
  cur <- read.table(fname, sep = ",", header = TRUE)
  if (s == min(subs)) Clickdata <- cur else Clickdata <- rbind(Clickdata, cur)
}
rm(cur)
Clickdata <- Clickdata %>% 
  filter(practiceNr == 3) %>% 
  mutate(uniqueID = paste0(subjectNr, "_", trialNr)) %>% select(-practiceNr) %>% 
  mutate(ymd = date_parse(str_extract(SysTime, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}"), format="%m/%d/%y"),
         SysTime = as_naive_time(year_month_day(get_year(ymd),get_month(ymd), get_day(ymd),
                                                SysTime_h, SysTime_m, SysTime_s, SysTime_ms,
                                                subsecond_precision = "millisecond"))) %>% 
  select(-c(SysTime_h, SysTime_m, SysTime_s, SysTime_ms, ymd))
Clickdata <- left_join(Clickdata,
                       ETdata %>% group_by(uniqueID) %>% slice(1) %>% 
                         rename(trialStartTime = SysTime) %>% select(c(trialStartTime, uniqueID))) %>% 
  group_by(uniqueID) %>% mutate(t = SysTime - trialStartTime) %>% ungroup() %>% 
  select(order(colnames(.)))

# =========================
# Plane projection (como en Ben: carga el preprocesado)
# =========================
load(raw_path("preProcessedETdata.Rdata"))  # añade planeIntersect_*

# =========================
# Exclusiones (igual que Ben)
# =========================
ETdata <- ETdata %>% group_by(uniqueID) %>% 
  mutate(invalidData = ifelse(planeIntersect_x < -1.5 | planeIntersect_x > 1.5 |
                                planeIntersect_y < -1.5 | planeIntersect_y > 1.5, 1, invalidData)) %>% 
  ungroup() %>% mutate(invalidData = replace_na(invalidData,1)) %>% 
  group_by(uniqueID) %>% 
  mutate(invalidWindow = slide_index_int(invalidData, t, max, .before=30, .after=30),
         planeIntersect_x = ifelse(invalidWindow==1, NA, planeIntersect_x),
         planeIntersect_y = ifelse(invalidWindow==1, NA, planeIntersect_y)) %>% 
  ungroup() %>% select(order(colnames(.)))

tmp <- ETdata %>% group_by(subjectNr, trialNr, uniqueID) %>% summarize(propInvalid = sum(invalidWindow)/n())
thresh <- .1
exclude_trials <- tmp %>% ungroup() %>% filter(propInvalid > thresh) %>% .$uniqueID
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)

tmp <- ETdata %>% group_by(subjectNr) %>% summarise(propInvalidTrials = 1-length(unique(trialNr))/24)
thresh <- .2
exclude_subs <- tmp %>% ungroup() %>% filter(propInvalidTrials > thresh) %>% .$subjectNr
ETdata <- ETdata %>% filter(!(subjectNr %in% exclude_subs))

totalTimes <- ETdata %>% group_by(uniqueID) %>% summarise(totalTime = max(t))
median_time <- median(totalTimes$totalTime)
longCutoff <- 2*median_time; shortCutoff <- median_time/2
long_trials  <- totalTimes$uniqueID[which(totalTimes$totalTime > longCutoff)] 
short_trials <- totalTimes$uniqueID[which(totalTimes$totalTime < shortCutoff)] 
exclude_trials <- c(short_trials, long_trials)
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)

tmp <- ETdata %>% group_by(uniqueID) %>% slice(1) %>% select(score, uniqueID)
scoreCutoff <- 6
exclude_trials <- tmp %>% filter(score <= scoreCutoff) %>% .$uniqueID
ETdata <- ETdata %>% filter(!(uniqueID %in% exclude_trials) | beepDetectionTrial==1)

# =========================
# Interpolación (igual que Ben) + **FIX de paréntesis**
# =========================
if (exists("Interpdata")) rm(Interpdata)
for (id in unique(ETdata$uniqueID)) {
  cur <- ETdata %>% filter(uniqueID == id)
  interp_points <- seq(0, max(cur$t), 1000/120)
  planeIntersect_x <- pracma::interp1(x = cur$t, y = cur$planeIntersect_x, xi = interp_points, method="linear")
  planeIntersect_y <- pracma::interp1(x = cur$t, y = cur$planeIntersect_y, xi = interp_points, method="linear")
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
  ungroup() %>% select(order(colnames(.)))

Interpdata <- Interpdata %>% 
  group_by(uniqueID) %>% 
  mutate(xdiff = c(0, diff(planeIntersect_x)),
         ydiff = c(0, diff(planeIntersect_y)),
         eucDist = sqrt(xdiff^2 + ydiff^2),
         euc_vel = eucDist/dt,
         theta2D = rad2deg(atan2(ydiff, xdiff)),
         ang_vel = theta2D/dt) %>%     # <<<<<< FIX: un solo paréntesis aquí
  select(-c(xdiff, ydiff, eucDist, theta2D)) %>% 
  ungroup() %>% 
  select(order(colnames(.)))

# =========================
# Guardado (mismos nombres que Ben, en data/processed)
# =========================
save(Trialdata, file=file.path(PROC_DIR, "Trialdata.Rdata")); write.csv(Trialdata, file.path(PROC_DIR,"Trialdata.csv"), row.names=FALSE)
save(Clickdata,  file=file.path(PROC_DIR, "Clickdata.Rdata")); write.csv(Clickdata,  file.path(PROC_DIR,"Clickdata.csv"),  row.names=FALSE)
save(ETdata,     file=file.path(PROC_DIR, "ETdata.Rdata"));    write.csv(ETdata,     file.path(PROC_DIR,"ETdata.csv"),     row.names=FALSE)
save(Interpdata, file=file.path(PROC_DIR, "Interpdata.Rdata"));write.csv(Interpdata, file.path(PROC_DIR,"Interpdata.csv"), row.names=FALSE)
write.csv(exclusionTally, file=file.path(PROC_DIR, "exclusionTally.csv"), row.names=FALSE)

message("✅ Finished. Inputs: data/raw ; Outputs: data/processed")
