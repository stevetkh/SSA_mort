library(rdhs)
library(demogsurv)
library(survival)
library(dplyr)

#' ## Steps in script
#'
#' 1. Extract sibling dataset from DHS
#' 2. Extract sibling dataset from MICS
#' 3. Calculate event & PYs data set 

#' ## DHS data sets

countries <- dhs_countries() %>%
  filter(RegionName == "Sub-Saharan Africa")
surveys <- dhs_surveys(countryIds = countries$DHS_CountryCode, surveyCharacteristicIds = 1)

#' Erroneously coded as SurveyCharacteristicId == 1 (maternal mortality)
#' * Ghana 1993 DHS, Sudan 1990 DHS
surveys <- filter(surveys, !SurveyId %in% c("GH1993DHS", "SD1990DHS"))

#' Note: the following surveys have sibling history data but datasets are restricted:
#' * Mauritania 2000-01 DHS
#' * Eritrea 1995 DHS


surveys <- surveys %>%
  left_join(
    select(countries, DHS_CountryCode, iso3 = ISO3_CountryCode, subregion = SubregionName),
    by = "DHS_CountryCode"
  ) %>%
  mutate(
    survey_id = paste0(iso3, SurveyYear, SurveyType),
    survey_year = as.integer(SurveyYear)
  )

ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", force = TRUE) %>%
  filter(SurveyId %in% surveys$SurveyId)
  
ird$path <- unlist(get_datasets(ird))

extract_dhs_sib <- function(path, SurveyId) {

  print(SurveyId)
  
  ir <- readRDS(path)
  sib <- reshape_sib_data(ir, widevars = grep("^v0|^aidsex", names(ir), value = TRUE))
  sib$SurveyId <- SurveyId
  
  if (substr(SurveyId, 1, 2) == "ET"){
    
    ## Convert Ethiopian calendar dates to Gregorian: add 92 months to CMC
    ## v007     "year of interview"
    ## v008     "date of interview (cmc)"
    ## v011     "date of birth (cmc)"
    ## mm4_01   "sibling's date of birth (cmc)"
    ## mm8_01   "date of death of sibling (cmc)"
    
    sib[c("v008", "v011", "mm4", "mm8")]  <- sib[c("v008", "v011", "mm4", "mm8")] + 92
    sib["v007"] <- 1900 + floor((sib["v008"] - 1) / 12)
  }

  if (inherits(sib$v007, "haven_labelled")) {
    print(paste0(SurveyId, ": v007 has labels:"))
    print(attr(sib$v007, "labels"))
  
    sib$v007 <- type.convert(haven::as_factor(sib$v007), as.is = TRUE)
  }

  sibrecode <- sib %>%
    transmute(
      SurveyId,
      stratum = paste(haven::as_factor(sib$v024), haven::as_factor(sib$v025)),
      cluster = v001,
      id = caseid,
      intvcmc = v008,
      intvy = v007,
      respsex = "female",
      respdobcmc = v011,
      weight = v005 / 1e6,
      sibidx = mmidx,
      sibsex = recode(as.integer(mm1), `1` = "male", `2` = "female", .default = NA_character_),
      sibdob = mm4,
      sibdod = mm8,
      sibdeath = recode(as.integer(sib$mm2), `0` = 1L, `1` = 0L, .default = NA_integer_)
    )

  sibrecode
}

sib_dhs_list <- Map(extract_dhs_sib, ird$path, ird$SurveyId)

sib_dhs <- sib_dhs_list %>%
  bind_rows() %>%
  left_join(
    surveys %>%
      select(iso3, country = CountryName, survey_id, subregion, SurveyId, survey_year, survey_type = SurveyType),
    by = "SurveyId"
  ) %>%
  select(iso3, country, survey_id, subregion, SurveyId, survey_year, survey_type, everything())

#' DHS male recode datasets with sibling history
#' * Zimbabwe 2005-06
#' * Zambia 2007
#'
#' Add "MR" to the survey_id

surveys_mr <- dhs_surveys(surveyIds = c("ZM2007DHS", "ZW2005DHS")) %>%
  left_join(
    select(countries, DHS_CountryCode, iso3 = ISO3_CountryCode, subregion = SubregionName),
    by = "DHS_CountryCode"
  ) %>%
  mutate(
    survey_id = paste0(iso3, SurveyYear, SurveyType, "MR"),
    survey_year = as.integer(SurveyYear)
  )

mrd <- dhs_datasets(fileType = "MR", fileFormat = "flat", force = TRUE) %>%
  filter(SurveyId %in% surveys_mr$SurveyId)
  
mrd$path <- unlist(get_datasets(mrd))

extract_dhs_sib_mr <- function(path, SurveyId) {

  print(SurveyId)
  
  mr <- readRDS(path)
  sib <- reshape_sib_data(mr, widevars = grep("^mv0", names(mr), value = TRUE),
                          idvar = "mcaseid", sib_idvar = "smidx", sibvar_regex = "^smm[idx0-9]")
  sib <- filter(sib, !is.na(smm1))
  sib$SurveyId <- SurveyId


  if (SurveyId == "ZW2005DHS") {

    ## ZW2005DHS has sibling age / age at death / years ago died only
    ## Note: This does NOT handle imputation for cases with missing data for
    ##   sibling current age, sibling age at death, or sibling died years ago.
    ##   --> missing values will result in slight underestimation of mortality
    ##  rates.
    
    sib <- sib %>%
      mutate(
        smm4 = case_when(smm2 == 1 & smm3 < 95 ~ mv008 - 12 * smm3 - 6,  # sibling alive
                        smm2 == 0 & smm6 < 95 & smm7 < 95 ~ mv008 - 12 * (smm6 + smm7) - 11), # sibling dead
        smm8 = case_when(smm2 == 0 & smm6 < 95 ~ mv008 - 12 * smm6 - 6)
      )
  }
  
  sibrecode <- sib %>%
    transmute(
      SurveyId,
      stratum = paste(haven::as_factor(sib$mv024), haven::as_factor(sib$mv025)),
      cluster = mv001,
      id = mcaseid,
      intvcmc = mv008,
      intvy = mv007,
      respsex = "male",
      respdobcmc = mv011,
      weight = mv005 / 1e6,
      sibidx = smidx,
      sibsex = recode(as.integer(smm1), `1` = "male", `2` = "female", .default = NA_character_),
      sibdob = smm4,
      sibdod = smm8,
      sibdeath = recode(as.integer(sib$smm2), `0` = 1L, `1` = 0L, .default = NA_integer_)
    )

  sibrecode
}

sib_dhs_mr_list <- Map(extract_dhs_sib_mr, mrd$path, mrd$SurveyId)

sib_dhs_mr <- sib_dhs_mr_list %>%
  bind_rows() %>%
  left_join(
    surveys_mr %>%
      select(iso3, country = CountryName, survey_id, subregion, SurveyId, survey_year, survey_type = SurveyType),
    by = "SurveyId"
  ) %>%
  select(iso3, country, survey_id, subregion, SurveyId, survey_year, survey_type, everything())

#' ## Add MICS surveys
#'
#' The following MICS surveys have maternal mortality module

## ben2014mics
## cog2014mics
## gin2016mics
## gnb2014mics
## mdg2018mics
## mrt2011mics
## mwi2013mics
## stp2014mics
## zwe2014mics
## zwe2019mics

surveys_mics_sib <- tribble(
  ~iso3, ~survey_year, ~sibsexvar, ~sibdobvar, ~sibdodvar, ~sibalivevar,
  "BEN",        2014,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "COG",        2014,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "GIN",        2016,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "GNB",        2014,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "MDG",        2018,     "MM15",    "MM17C",   "MM18C",        "MM16",
  "MRT",        2011,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "MWI",        2013,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "STP",        2014,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "ZWE",        2014,      "MM5",     "MM7C",    "MM8C",         "MM6",
  "ZWE",        2019,     "MM15",    "MM17C",   "MM18C",        "MM16") %>%
  left_join(
    select(countries, iso3 = ISO3_CountryCode, country = CountryName, subregion = SubregionName) %>%
      bind_rows(data.frame(iso3 = "GNB", country = "Guinea-Bissau", subregion = "Western Africa"))
  ) %>%
  mutate(
    survey_type = "MICS",
    survey_id = paste0(iso3, survey_year, survey_type),
    path = paste0("raw/mics/", tolower(survey_id), ".rds")
  )
  
  
#' Standard MICS variable codings:
#' - They are different for the 2 surveys conducted after 2018
#' - A couple of surveys have WMWEIGHT instead of wmweight
#' - One survey is missing WDOB variable; needs to be merged from
#'   women's recode if it is needed.
#' 
#' * HH1      "Cluster number"
#' * HH2      "Household number"
#' * LN       "Line number"
#' * MMLN     "Sibling's line number"
#' * MM5      "Sibling's gender"
#' * MM6      "Sibling still alive"
#' * MM7      "Age of sibling"
#' * MM8      "Years since death"
#' * MM9      "Age at death of sibling"
#' * MM10     "Pregnant when died"
#' * MM11     "Died during the childbirth"
#' * MM12     "Died within two months after the end of a pregnancy or childbirth"
#' * MM13     "Number of live born children during lifetime"
#' * HH6      "Area"
#' * HH7      "Province"
#' * WDOI     "Date of interview women (CMC)"
#' * WDOB     "Date of birth of woman (CMC)"
#' * MM7C     "Imputed date of birth"
#' * MM8C     "Imputed date of death"
#' * welevel  "Education"
#' * wmweight "Women's sample weight"
#' * wscore   "Combined wealth score"
#' * windex5  "Wealth index quintile"
#' * wscoreu  "Urban wealth score"
#' * windex5u "Urban wealth index quintile"
#' * wscorer  "Rural wealth score"
#' * windex5r "Rural wealth index quintile"


extract_mics_sib <- function(path, survey_id, strata=c("HH7", "HH6"),
                             sibsexvar = "MM5", sibdobvar = "MM7C", sibdodvar = "MM8C",
                             sibalivevar = "MM6") {

  print(survey_id)
  
  mcs <- readRDS(path)
  sib <- mcs$mm

  if ( ! "WDOB" %in% names(sib)) {
    ## For MDG2018MICS which is missing WDOB in mm dataset
    sib$WDOB <- NA_integer_
  }

  if ("WMWEIGHT" %in% names(sib)) {
    ## For MDG2018MICS which is missing WDOB in mm dataset
    sib$wmweight <- sib$WMWEIGHT
  }

  sibrecode <- sib %>%
    transmute(
      survey_id = survey_id,
      stratum = do.call(paste, haven::as_factor(sib[strata])),
      cluster = HH1,
      id = paste(HH1, HH2, LN),
      intvcmc = WDOI,
      intvy = floor((WDOI - 1)/12 + 1900),
      respsex = "female",
      respdobcmc = WDOB,
      weight = wmweight,
      sibidx = MMLN,
      sibsex = recode(as.integer(sib[[sibsexvar]]),
                      `1` = "male", `2` = "female", .default = NA_character_),
      sibdob = sib[[sibdobvar]],
      sibdod = sib[[sibdodvar]],
      sibdeath = recode(as.integer(sib[[sibalivevar]]),
                        `2` = 1L, `1` = 0L, .default = NA_integer_)
    )

  sibrecode
}

sib_mics_list <- Map(extract_mics_sib,
                     path = surveys_mics_sib$path,
                     survey_id = surveys_mics_sib$survey_id,
                     sibsexvar = surveys_mics_sib$sibsexvar,
                     sibdobvar = surveys_mics_sib$sibdobvar,
                     sibdodvar = surveys_mics_sib$sibdodvar,
                     sibalivevar = surveys_mics_sib$sibalivevar)

sib_mics <- sib_mics_list %>%
  bind_rows() %>%
  left_join(
    surveys_mics_sib %>%
      select(iso3, country, survey_id, subregion, survey_year, survey_type),
    by = "survey_id"
  ) %>%
  select(iso3, country, survey_id, subregion, survey_year, survey_type, everything())

#' ## Pool all sibling data and save

sib_data <- bind_rows(sib_dhs, sib_dhs_mr, sib_mics)

#' Remove 2 rows with missing weight from STP 2014 MICS
sib_data %>%
  filter(is.na(weight))

sib_data <- filter(sib_data, !is.na(weight))


#' Recode as factors, which demogsurv requires for stratification variables

sib_data <- sib_data %>%
  mutate(
    iso3 = factor(iso3),
    country = factor(country),
    survey_id = factor(survey_id),
    subregion = factor(subregion),
    survey_type = factor(survey_type),
    stratum = factor(stratum),
    respsex = factor(respsex),
    sibsex = factor(sibsex)
  )
    
dir.create("data")
saveRDS(sib_data, "data/sib_data.rds")

sib_surveys <- distinct(sib_data, iso3, country, survey_id, subregion, SurveyId, survey_year, survey_type)
saveRDS(sib_surveys, "data/sib_surveys.rds")

#' ## Deaths and person-years dataset 

sib_data_spl <- sib_data %>%
  mutate(tstop = if_else(sibdeath == 1, sibdod, intvcmc)) %>%
  split(~survey_id)

pyears_data_spl <- lapply(sib_data_spl, 
                          function(x){
                            print(as.character(x$survey_id[1]))
                            demog_pyears(~survey_id + sibsex,
                                         data = x,
                                         period = -15:1 + x$survey_year[1],
                                         agegr = 0:70,
                                         tips = 0:15,
                                         dob = "sibdob",
                                         intv = "intvcmc",
                                         event = "sibdeath",
                                         tstart = "sibdob",
                                         tstop = "tstop",
                                         weights = "weight")$data
                          })

pyears_data <- pyears_data_spl %>%
  bind_rows() %>%
  right_join(sib_surveys, ., by = "survey_id") %>%
  select(names(sib_surveys), everything())

saveRDS(pyears_data, "data/pyears_data.rds")


#' 
