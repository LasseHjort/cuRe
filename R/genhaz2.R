general.haz2 <- function (formula, data, weights, subset, na.action, rmap,
                          ratetable = survexp.dk, scale = 1, model = FALSE,
                          x = FALSE, y = FALSE){
  Call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m <- m[c(1, match(c("formula", "data", "weights", "subset",
                      "na.action"), names(m), nomatch = 0))]
  m[[1L]] <- quote(stats::model.frame)
  Terms <- if (missing(data))
    terms(formula, "ratetable")
  else terms(formula, "ratetable", data = data)
  rate <- attr(Terms, "specials")$ratetable
  if (length(rate) > 1)
    stop("Can have only 1 ratetable() call in a formula")
  if (length(rate) == 1) {
    if (!missing(rmap))
      stop("The ratetable() call in a formula is depreciated")
    stemp <- untangle.specials(Terms, "ratetable")
    rcall <- as.call(parse(text = stemp$var)[[1]])
    rcall[[1]] <- as.name("list")
    Terms <- Terms[-stemp$terms]
  }
  else if (!missing(rmap)) {
    rcall <- substitute(rmap)
    if (!is.call(rcall) || rcall[[1]] != as.name("list"))
      stop("Invalid rcall argument")
  }
  else rcall <- NULL
  if (is.ratetable(ratetable))
    varlist <- attr(ratetable, "dimid")
  else if (inherits(ratetable, "coxph")) {
    varlist <- all.vars(delete.response(ratetable$terms))
  }
  else stop("Invalid rate table")
  temp <- match(names(rcall)[-1], varlist)
  if (any(is.na(temp)))
    stop("Variable not found in the ratetable:", (names(rcall))[is.na(temp)])
  if (any(!(varlist %in% names(rcall)))) {
    to.add <- varlist[!(varlist %in% names(rcall))]
    temp1 <- paste(text = paste(to.add, to.add, sep = "="),
                   collapse = ",")
    if (is.null(rcall))
      rcall <- parse(text = paste("list(", temp1, ")"))[[1]]
    else {
      temp2 <- deparse(rcall)
      rcall <- parse(text = paste("c(", temp2, ",list(",
                                  temp1, "))"))[[1]]
    }
  }
  newvar <- all.vars(rcall)
  if (length(newvar) > 0) {
    tform <- paste(paste(deparse(Terms), collapse = ""),
                   paste(newvar, collapse = "+"), sep = "+")
    m$formula <- as.formula(tform, environment(Terms))
  }
  m <- eval(m, parent.frame())
  n <- nrow(m)
  if (n == 0)
    stop("Data set has 0 rows")
  weights <- model.extract(m, "weights")
  if (length(weights) == 0)
    weights <- rep(1, n)
  if (class(ratetable) == "ratetable" && any(weights != 1))
    warning("weights ignored")
  if (any(attr(Terms, "order") > 1))
    stop("Survexp cannot have interaction terms")
  Y <- model.extract(m, "response")
  if (is.matrix(Y)) {
    if (is.Surv(Y) && attr(Y, "type") == "right")
      Y <- Y[, 1]
    else stop("Illegal response value")
  }
  if (any(Y < 0))
    stop("Negative follow up time")
  newtime <- Y
  ovars <- attr(Terms, "term.labels")
  rdata <- data.frame(eval(rcall, m), stringsAsFactors = TRUE)
  if(!is.ratetable(ratetable)) stop("Invalid ratetable")
  rtemp <- match.ratetable(rdata, ratetable)
  R <- rtemp$R
  if (TRUE) {
    if (israte)
      temp <- survexp.fit(1:n, R, Y, max(Y), TRUE, ratetable)
    else {
      rmatch <- match(names(data), names(rdata))
      if (any(is.na(rmatch)))
        rdata <- cbind(rdata, data[, is.na(rmatch)])
      temp <- survexp.cfit(1:n, rdata, Y, "individual",
                           ratetable)
    }
    if (method == "individual.s")
      xx <- temp$surv
    else xx <- -log(temp$surv)
    names(xx) <- row.names(m)
    na.action <- attr(m, "na.action")
    if (length(na.action))
      return(naresid(na.action, xx))
    else return(xx)
  }
  if (length(ovars) == 0)
    X <- rep(1, n)
  else {
    odim <- length(ovars)
    for (i in 1:odim) {
      temp <- m[[ovars[i]]]
      ctemp <- class(temp)
      if (!is.null(ctemp) && ctemp == "tcut")
        stop("Can't use tcut variables in expected survival")
    }
    X <- strata(m[ovars])
  }
  if (israte)
    temp <- survexp.fit(as.numeric(X), R, Y, newtime, method ==
                          "conditional", ratetable)
  else {
    temp <- survexp.cfit(as.numeric(X), rdata, Y, method,
                         ratetable, weights)
    newtime <- temp$time
  }
  if (missing(times)) {
    n.risk <- temp$n
    surv <- temp$surv
  }
  else {
    if (israte)
      keep <- match(times, newtime)
    else {
      n <- length(temp$time)
      keep <- approx(temp$time, 1:n, xout = times, yleft = 0,
                     method = "constant", f = 0, rule = 2)$y
    }
    if (is.matrix(temp$surv)) {
      surv <- (rbind(1, temp$surv))[keep + 1, , drop = FALSE]
      n.risk <- temp$n[pmax(1, keep), , drop = FALSE]
    }
    else {
      surv <- (c(1, temp$surv))[keep + 1]
      n.risk <- temp$n[pmax(1, keep)]
    }
    newtime <- times
  }
  newtime <- newtime/scale
  if (is.matrix(surv)) {
    dimnames(surv) <- list(NULL, levels(X))
    out <- list(call = Call, surv = drop(surv), n.risk = drop(n.risk),
                time = newtime)
  }
  else {
    out <- list(call = Call, surv = c(surv), n.risk = c(n.risk),
                time = newtime)
  }
  if (model)
    out$model <- m
  else {
    if (x)
      out$x <- X
    if (y)
      out$y <- Y
  }
  if (israte && !is.null(rtemp$summ))
    out$summ <- rtemp$summ
  if (no.Y)
    out$method <- "Ederer"
  else if (conditional)
    out$method <- "conditional"
  else out$method <- "cohort"
  class(out) <- c("survexp", "survfit")
  out
}
