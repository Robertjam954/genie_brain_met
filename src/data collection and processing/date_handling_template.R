############################################################
# 📘 DATE HANDLING TEMPLATE — Working With Dates and Times #
# Section 3.6 — Lubridate & nycflights13 Example
############################################################

# -----------------------------
# 1. Setup
# -----------------------------
# Install (if needed) and load packages
# install.packages('lubridate')
# install.packages('nycflights13')
library(lubridate)
library(dplyr)
library(nycflights13)

# -----------------------------
# 2. Creating Date and Date-Time Objects
# -----------------------------
# --- From strings ---

# Year-Month-Day format
date_ymd <- ymd("1988-09-29")

# Month-Day-Year format
date_mdy <- mdy("September 29th, 1988")

# Day-Month-Year format
date_dmy <- dmy("29-Sep-1988")

# Display results
date_ymd; date_mdy; date_dmy

# --- From strings including time ---
datetime_ymd_hms <- ymd_hms("1988-09-29 20:11:59")
datetime_ymd_hms

# -----------------------------
# 3. Creating Dates/Date-Times from Individual Parts
# -----------------------------

# --- Using make_date() ---
flights_dates <- flights %>%
  select(year, month, day) %>%
  mutate(departure = make_date(year, month, day))

head(flights_dates)

# --- Using make_datetime() ---
flights_datetimes <- flights %>%
  select(year, month, day, hour, minute) %>%
  mutate(departure = make_datetime(year, month, day, hour, minute))

head(flights_datetimes)

# -----------------------------
# 4. Extracting Components from Dates
# -----------------------------

mydate <- ymd("1988-09-29")

# Extract year
year_info <- year(mydate)

# Extract day of the month
day_info <- mday(mydate)

# Extract weekday (numeric)
weekday_num <- wday(mydate)

# Extract weekday (labeled)
weekday_label <- wday(mydate, label = TRUE)

# Display results
year_info; day_info; weekday_num; weekday_label

# -----------------------------
# 5. Working with Time Spans
# -----------------------------

# Example: calculate age
birthday <- ymd("1988-09-29")

# Subtract from today's date
age_diff <- today() - birthday
age_diff   # in days

# Convert to duration (in seconds, then years)
age_duration <- as.duration(age_diff)
age_duration

# Approximate years old
years_old <- as.numeric(age_duration) / (60 * 60 * 24 * 365)
years_old

# -----------------------------
# 6. Additional Examples
# -----------------------------
# Adding or subtracting days/months/years
today_plus_30 <- today() + days(30)
today_minus_1year <- today() - years(1)

today_plus_30
today_minus_1year

############################################################
# ✅ END OF TEMPLATE
# Use this file as a quick reference or teaching resource
############################################################
