#' Extract general population hazard
#'
#' Yearly general population hazards matched on age, gender, and calendar year is extracted from a ratetable.
#'
#' @param time Either a numeric vector of follow-up times (in days) or a character indicating the variable
#' containing the follow-up times in the data.
#' @param rmap A named list. The names must match the dimension names of the ratetable input.
#' The elements should be given as vectors containing the relevant variables in the data
#' or a character indicating the name of the variable in the data.
#' Elements in the list named 'age' and 'year' are transformed such that the age and year of
#' the event/censoring time is used.
#' 'age' must be a numeric vector given as age in days
#' or a character indicating the variable containing the patient ages in the data.
#' 'year' must be of class \code{Date} with the calendar time points
#' or a character indicating the variable containing the calendar times in the data.
#' Other variables should be coded as characters in the data.
#' @param data The data from which to extract variables from.
#' If \code{time}, \code{age}, \code{sex}, or \code{year} are not characters, this will not be used.
#' @param ratetable Object of class \code{ratetable} to extract background hazards from. Defaults to \code{survexp.dk}.
#' @param scale Numeric to adjust the scale of the outputted hazard values.
#' If the ratetable provides daily hazards and \code{scale = 365.24} (default), the outputted hazard values are yearly.
#' @return An object of class \code{numeric} containing the yearly expected hazards.
#' @export
#' @example inst/general.haz.ex.R
#' @importFrom reshape2 melt


general.haz <- function(time, rmap, data = NULL, ratetable = cuRe::survexp.dk, scale = 365.24){

  #Compute order of the dimnames in the ratetable - may not match with the attributes
  dimid <- attr(ratetable, "dimid")
  if(is.null(dimid)){
    dimid <- names(attr(ratetable, "dimnames"))
  }
  od <- sapply(names(rmap), function(x) which(dimid == x))
  n <- length(time)

  #Stop if not all dimension names in "ratetable" input are given in "rmap" input
  for(var in dimid){
    if(!(var %in% names(rmap))) stop('Not all dimension names in ratetable specified in "rmap" argument')
  }

  #Collect "rmap" covariate values from "data" in one data frame.
  #Extract those values from data if only characters are provided.
  if(is.character(time)){
    time <- data[, time]
  }

  rmap_vars <- as.data.frame(matrix(nrow=length(time), ncol=length(rmap)))
  names(rmap_vars) <- names(rmap)

  for(varname in names(rmap)) {
    if(length(rmap[varname][[1]])==1){ #True if input given as character
      rmap_vars[, varname] <- data[, rmap[varname][[1]]]
    } else {
      rmap_vars[, varname] <- rmap[varname][[1]]
    }
    if(is.na(match(varname, c("age", "year")))){
      rmap_vars[, varname] <- as.character(rmap_vars[, varname])
    }
  }

  #Transformation of age to age at event/censoring time
  if("age" %in% names(rmap)){
    #Compute max age in the ratetable
    max_age <- max(as.numeric(attr(ratetable, "dimnames")[[od["age"]]]))

    #Compute age and year after "time" days
    age_new <- pmin(round((rmap_vars[, "age"] + time) / ayear), max_age)
    rmap_vars[, "age"] <- age_new
  }

  #If years expand that available in the ratetable take the closest available data points
  if("year" %in% names(rmap)){
    year_eval <- format(rmap_vars[, "year"] + time, "%Y")

    ryear <- range(as.numeric(dimnames(ratetable)[[od["year"]]]))
    year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
    year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)

    rmap_vars[, "year"] <- year_eval
  }

  #Create useful data frame from "ratetable" input and use dimnames from the ratetable
  df <- melt(as.matrix(ratetable), as.is = T)
  names(df)[1:length(dimid)] <- dimid

  #Order the dataframe to match the dimensions from the ratetable
  rmap_vars <- rmap_vars[, od]
  #Define order of individuals to go back after merging
  rmap_vars$ord <- 1:nrow(rmap_vars)

  #Merge the individual data frame with melted ratetable - only keep rows in individual data frame
  haz_df <- merge(rmap_vars, df, by = names(rmap), all.x = T)

  #Extract ratetable values and use same order as input
  haz <- haz_df$value[order(haz_df$ord)]

  #Multiple by scale and output
  haz * scale
}

# attr <- attributes(ratetable)
# attr$dim
#
# a <- survexp.dk[,,"male"]
# b <- survexp.dk[,,"female"]
# a <- melt(as.matrix(a))
#
#
# b <- melt(as.matrix(b))
# a$sex <- "male"
# b$sex <- "female"
#
# attr$dimid
# df <- rbind(a, b)
#
# names(df)[1:2] <- attr$dimid[1:2]
#
#
# d <- merge(df, D, by = c("age", "sex", "year"), all.y = T)
# d$value


# general.haz2 <- function(time, rmap, data = NULL, ratetable = survexp.dk, scale = ayear){
#   dimid <- attr(ratetable, "dimid")
#   vars <-
#
#     od <- sapply(c("age", "sex", "year"), function(x) which(dimid == x))
#   n <- length(time)
#
#   haz <- rep(NA, n)
#   sex_new <- as.character(sex)
#   age_new <- pmin(round((age + time) / ayear), 99)
#   year_eval <- format(year + time, "%Y")
#   ryear <- range(as.numeric(dimnames(ratetable)[[od["year"]]]))
#   year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
#   year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)
#
#
#   D <- data.frame(age = age_new, sex = sex_new, year = year_eval, stringsAsFactors = F)
#   D <- D[, od]
#
#   #D2 <- D
#   #D2$sex <- ifelse(D2$sex == "male", 1, 2)
#   #D2$year <- as.numeric(D2$year)
#   #a <- match.ratetable(as.matrix(D2), ratetable)
#
#   dim_names <- dimnames(ratetable)
#   J <- data.frame(rates = c(ratetable), rep(as.numeric(dim_names[[1]]), length(dim_names[[2]]) * length(dim_names[[3]])),
#                   rep(as.numeric(dim_names[[2]]), length(dim_names[[1]]) * length(dim_names[[3]])),
#                   rep(dim_names[[3]], each = length(dim_names[[1]]) *  length(dim_names[[2]])))
#
#   names(J)[-1] <- dimid
#
#   #head(J)
#
#   #ratetable["38", "1835",]
#
#   merge(D, J)$rates * scale
# }
