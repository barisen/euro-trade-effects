# =================================================================
# EURO EFFECT - ROSE EFFECT
# CEPII GRAVITY MODAL
# Data can be downloaded from: https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8
# =================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})
setFixest_nthreads(1)

# --------------------------------------------------------
# 1. CONFIGURATION
# --------------------------------------------------------
DATA_PATH <- "/Users/barissen/Documents/Personal/HU_Master/Foreign_Trade_25/Final_Paper/Data/Gravity_csv_V202211/Gravity_V202211.csv"
FLOW_COL  <- "tradeflow_imf_o"
SAMPLE       <- "EMU_plusROW"
ROW_FILTER   <- "allROW"
Y_START <- 1980
Y_END   <- 2021
MIN_FLOW <- 5000
EARLY_EMU_ONLY <- TRUE
INCLUDE_DE <- TRUE
INCLUDE_NL <- TRUE
INCLUDE_GRC <- TRUE
INCLUDE_ITA <- TRUE

# --------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# --------------------------------------------------------
message("Loading data from: ", DATA_PATH)
if (!file.exists(DATA_PATH)) stop("Data file not found.")
dt <- fread(DATA_PATH)

if (!("country_id_o" %in% names(dt))) dt[, country_id_o := as.integer(factor(iso3_o))]
if (!("country_id_d" %in% names(dt))) dt[, country_id_d := as.integer(factor(iso3_d))]

dt[, trade := as.numeric(get(FLOW_COL))]
dt <- dt[is.finite(trade) & trade >= MIN_FLOW]

# --------------------------------------------------------
# 3. VARIABLE CREATION
# --------------------------------------------------------
euro_year_map <- c(AUT=1999, BEL=1999, FIN=1999, FRA=1999, DEU=1999, IRL=1999,
                   ITA=1999, LUX=1999, NLD=1999, PRT=1999, ESP=1999, GRC=2001)

dt[, euro_year_o := euro_year_map[iso3_o]]
dt[, euro_year_d := euro_year_map[iso3_d]]
dt[is.na(euro_year_o), euro_year_o := 9999L]
dt[is.na(euro_year_d), euro_year_d := 9999L]

dt[, euro_pair := as.integer(year >= euro_year_o & year >= euro_year_d & euro_year_o < 9999 & euro_year_d < 9999)]
dt[, pre_euro_pair := as.integer((year %in% 1996:1998) & (euro_year_o <= 2001 & euro_year_d <= 2001))]
dt[, treatment_year := 0L]
dt[euro_year_o <= 2001 & euro_year_d <= 2001, treatment_year := 1999L]

# --------------------------------------------------------
# 4. FILTERING AND PREP
# --------------------------------------------------------
emu_set <- names(euro_year_map)
if (SAMPLE == "EMU_intra") {
  d <- dt[iso3_o %in% emu_set & iso3_d %in% emu_set & iso3_o != iso3_d]
} else {
  d <- dt[(iso3_o %in% emu_set | iso3_d %in% emu_set) & iso3_o != iso3_d]
}
d <- d[year >= Y_START & year <= Y_END]
d[, pair_id := paste(pmin(iso3_o, iso3_d), pmax(iso3_o, iso3_d), sep = "_")]

message(sprintf("Toggles -> DE:%s NL:%s GRC:%s ITA:%s | EARLY_EMU_ONLY:%s",
                INCLUDE_DE, INCLUDE_NL, INCLUDE_GRC, INCLUDE_ITA, EARLY_EMU_ONLY))
message("Rows in estimation sample: ", format(nrow(d), big.mark=","))

# ---- Helper functions ----
add_exporter_interactions <- function(DT, all_trade = FALSE, use_DE=FALSE, use_NL=FALSE, use_GRC=FALSE, use_ITA=FALSE) {
  side <- function(code) if (!all_trade) paste0("(iso3_o == \"", code, "\")") else paste0("((iso3_o == \"", code, "\") | (iso3_d == \"", code, "\"))")
  if (use_DE)  DT[, euro_DE  := euro_pair * eval(parse(text=side("DEU")))]
  if (use_NL)  DT[, euro_NL  := euro_pair * eval(parse(text=side("NLD")))]
  if (use_GRC) DT[, euro_GRC := euro_pair * eval(parse(text=side("GRC")))]
  if (use_ITA) DT[, euro_ITA := euro_pair * eval(parse(text=side("ITA")))]
  invisible(DT)
}
rhs_from_toggles <- function(base="euro_pair", use_DE=FALSE, use_NL=FALSE, use_GRC=FALSE, use_ITA=FALSE) {
  terms <- base
  if (use_DE)  terms <- paste(terms, "euro_DE",  sep=" + ")
  if (use_NL)  terms <- paste(terms, "euro_NL",  sep=" + ")
  if (use_GRC) terms <- paste(terms, "euro_GRC", sep=" + ")
  if (use_ITA) terms <- paste(terms, "euro_ITA", sep=" + ")
  terms
}

# ========================================================
# PART A: ORIGINAL 5-MODEL ANALYSIS
# ========================================================
message("\n--- Running Original 5-Model Specification ---")

# --- Model 1: Exporter-only PPML ---
d1 <- copy(d)
add_exporter_interactions(d1, all_trade = FALSE, use_DE=INCLUDE_DE, use_NL=INCLUDE_NL, use_GRC=INCLUDE_GRC, use_ITA=INCLUDE_ITA)
rhs1 <- rhs_from_toggles("euro_pair", use_DE=INCLUDE_DE, use_NL=INCLUDE_NL, use_GRC=INCLUDE_GRC, use_ITA=INCLUDE_ITA)
fml1 <- as.formula(paste("trade ~", rhs1, "| pair_id + country_id_o^year + country_id_d^year"))
mod_exp <- fepois(fml1, cluster = ~ pair_id, data = d1)

# --- Model 2: All-trade PPML ---
d2 <- copy(d)
add_exporter_interactions(d2, all_trade = TRUE, use_DE=INCLUDE_DE, use_NL=INCLUDE_NL, use_GRC=INCLUDE_GRC, use_ITA=INCLUDE_ITA)
rhs2 <- rhs_from_toggles("euro_pair", use_DE=INCLUDE_DE, use_NL=INCLUDE_NL, use_GRC=INCLUDE_GRC, use_ITA=INCLUDE_ITA)
fml2 <- as.formula(paste("trade ~", rhs2, "| pair_id + country_id_o^year + country_id_d^year"))
mod_all <- fepois(fml2, cluster = ~ pair_id, data = d2)

# --- Model 3: Trade Surplus OLS (ENHANCED WITH PAIR FE) ---
exports_dt <- d[, .(exports = trade, country = iso3_o, partner = iso3_d, year)]
imports_dt <- d[, .(imports = trade, country = iso3_d, partner = iso3_o, year)]
surplus <- merge(exports_dt, imports_dt, by = c("country", "partner", "year"), all = TRUE)
surplus[is.na(exports), exports := 0]
surplus[is.na(imports), imports := 0]
surplus[, trade_surplus := exports - imports]
surplus[, euro_year_country := euro_year_map[country]]
surplus[, euro_year_partner := euro_year_map[partner]]
surplus[!is.na(euro_year_country) & !is.na(euro_year_partner), euro_pair := as.integer(year >= euro_year_country & year >= euro_year_partner)]
surplus[is.na(euro_pair), euro_pair := 0]
# NEW: Create a pair ID for the surplus data
surplus[, pair_cp := paste(pmin(country, partner), pmax(country, partner), sep = "_")]
if (INCLUDE_DE) surplus[, euro_DE := euro_pair * (country == "DEU")]
if (INCLUDE_NL) surplus[, euro_NL := euro_pair * (country == "NLD")]
if (INCLUDE_GRC) surplus[, euro_GRC := euro_pair * (country == "GRC")]
if (INCLUDE_ITA) surplus[, euro_ITA := euro_pair * (country == "ITA")]
rhs3 <- rhs_from_toggles("euro_pair", use_DE=INCLUDE_DE, use_NL=INCLUDE_NL, use_GRC=INCLUDE_GRC, use_ITA=INCLUDE_ITA)
# NEW: Add pair_cp to the fixed effects formula
fml3 <- as.formula(paste("trade_surplus ~", rhs3, "| pair_cp + country^year + partner^year"))
mod_surplus <- feols(fml3, data = surplus)

# --- Model 4: OLS log(trade) ---
d_ols <- d1[trade > 0]
d_ols[, log_trade := log(trade)]
fml4 <- as.formula(paste("log_trade ~", rhs1, "| pair_id + country_id_o^year + country_id_d^year"))
mod_ols_log <- feols(fml4, cluster = ~ pair_id, data = d_ols)

# --- Model 5: OLS asinh(trade) ---
d_asinh <- d1
d_asinh[, asinh_trade := asinh(trade)]
fml5 <- as.formula(paste("asinh_trade ~", rhs1, "| pair_id + country_id_o^year + country_id_d^year"))
mod_ols_asinh <- feols(fml5, cluster = ~ pair_id, data = d_asinh)

# ========================================================
# PART B: NEW ROBUSTNESS CHECKS
# ========================================================
message("\n--- Running New Robustness Checks ---")

# --- Model 6: Placebo Test (PPML) ---
mod_placebo <- fepois(trade ~ euro_pair + pre_euro_pair | pair_id + country_id_o^year + country_id_d^year,
                      data = d, cluster = ~pair_id)

# --- Model 7: Trimmed Event Study ---
p99 <- quantile(d$trade, 0.99, na.rm = TRUE)
d_trimmed <- d[trade < p99]
message(paste0("\nTrimming event study data at the 99th percentile: ", format(p99, big.mark=",")))
message(paste0("Rows in trimmed sample for event study: ", format(nrow(d_trimmed), big.mark=",")))
mod_event_study_trimmed <- fepois(trade ~ sunab(treatment_year, year, ref.p = -1) | pair_id + country_id_o^year + country_id_d^year,
                                  data = d_trimmed, cluster = ~pair_id)

# ========================================================
# PART C: OUTPUTS
# ========================================================
message("\n\n--- TABLE 1: Results from Original 5-Model Specification ---")
etable(
  mod_exp, mod_all, mod_surplus, mod_ols_log, mod_ols_asinh,
  headers = c("1) PPML Exp-Only", "2) PPML All-Trade", "3) OLS Surplus (3-way FE)", "4) OLS log(trade)", "5) OLS asinh(trade)"),
  dict = c(euro_pair="Euro (Avg)", euro_DE="Euro x DE", euro_NL="Euro x NL", euro_GRC="Euro x GRC", euro_ITA="Euro x ITA"),
  tex = FALSE
)

message("\n\n--- TABLE 2: Results from New Robustness Checks ---")
etable(
  mod_placebo,
  headers = c("Placebo Test (PPML)"),
  dict = c(euro_pair = "Actual Euro Effect", pre_euro_pair = "Placebo Effect (Pre-Euro)"),
  notes = "All models include Pair + Exporter-Year + Importer-Year fixed effects.",
  tex = FALSE
)

message("\n--- Generating Clean Event Study Plot (from trimmed data) ---")
iplot(mod_event_study_trimmed,
      xlab = 'Years relative to Euro Adoption',
      main = 'Dynamic Effect of the Euro (Top 1% of Trade Flows Removed)')
