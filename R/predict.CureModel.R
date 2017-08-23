
#' Predict function for cure models
#'
#' This function is used to make predictions of the cure models.
#'
#' @param fit An object of type CureModel.
#' @param newdata A dataframe with new data to predict from.
#' @param times A numeric vector including the time points for which to make predictions.
#' @param type A character denoting the type of prediction. Possible values are \code{relsurv} (default), \code{curerate}, and \code{probcure} (see details)
#' @param ci Logical indicating whether to include confidence intervals. Default is \code{TRUE}
#' @return An object of class \code{matrix} including the predictions.
#' @details For type \code{relsurv}, relative survival predictions are made, while for \code{curerate}, \eqn{\pi} and for \code{probcure}\cr
#' the probability\cr
#' \deqn{\p{\text{prob}| T > t} = \frac{\pi}{\pi + (1 - \pi) S_u(t)}}-

predict.CureModel <- function(fit, newdata, times, type = "relsurv", ci = T){
  link_fun <- get_link(fit$link)
  if(type == "relsurv"){
    surv_fun <- get_surv2(fit$dist)
    nr.coefs <- cumsum(unlist(lapply(fit$coefs, length)))
    nr.coefs_indi <- unlist(lapply(fit$coefs, length))
    rs_fun <- function(times, vars, lps = list(1, 1, 1, 1)){
      pi_terms <- link_fun(lps[[1]] %*% vars[1:nr.coefs[1]]) + (1 - link_fun(lps[[1]] %*% vars[1:nr.coefs[1]]))
      surv_term <- surv_fun(times,
                            k1 = vars[(nr.coefs[1] + 1):nr.coefs[2]],
                            k2 = vars[ifelse(nr.coefs_indi[3] != 0, (nr.coefs[2] + 1):nr.coefs[3], NA)],
                            k3 = vars[ifelse(nr.coefs_indi[4] != 0, (nr.coefs[3] + 1):nr.coefs[4], NA)],
                            lp1 = lps[[2]], lp2 = lps[[3]], lp3 = lps[[4]])
      pi_terms * surv_term
    }
    if(is.null(newdata)){
      vars <- c(link_fun(fit$coefs[[1]]), unlist(fit$coefs[-1]))
      rs <- rs_fun(times, vars = vars)
      rss <- data.frame(times = times, RS = rs)
      if(ci){
        grads <- jacobian(rs_fun, x = vars, times = times)
        VAR <- apply(grads, MARGIN = 1, function(x) x %*% fit$cov %*% x)
        rss$VAR <- VAR
        rss$ci.upper <- rss$RS + qnorm(0.975) * sqrt(rss$VAR)
        rss$ci.lower <- rss$RS - qnorm(0.975) * sqrt(rss$VAR)
      }
      rss
    }else{
      tt <- terms(fit$formulas$formula.gamma)
      formula.2 <- formula(delete.response(tt))
      fit$formulas$formula.gamma <- formula.2
      design_matrices <- lapply(fit$formulas, function(x){
        if(length(x) > 0){
          return(model.matrix(x, data = newdata))
        }
      })


      #lps <- vector("list", length(design_matrices))
      #for(i in 1:length(lps)){
      #  if(length(design_matrices[[i]]) > 0){
      #    lps[[i]] <- design_matrices[[i]] %*% fit$coefs[[i]]
      #  }
      #}

      #lps <- do.call(cbind, lps)
      rss <- matrix(nrow = length(times), ncol = nrow(newdata))
      vars <- unlist(fit$coefs)
      for(i in 1:nrow(newdata)){
        rss[, i] <- rs_fun(times = times, vars = vars, lps = lapply(design_matrices, function(x) x[i,]))
      }
      colnames(rss) <- paste0("RS", 1:ncol(rss))
      rss <- cbind(times, rss)
      rss
    }
  }else if(type == "curerate"){
    if(is.null(newdata)){
      pi <- link_fun(fit$coefs[[1]])
      VAR <- diag(fit$cov)[1:length(pi)]
      ci.upper <- link_fun(fit$coefs[[1]] + qnorm(0.975) * sqrt(VAR))
      ci.lower <- link_fun(fit$coefs[[1]] - qnorm(0.975) * sqrt(VAR))
      data.frame(pi, ci.upper, ci.lower)
    }else{
      link_fun <- get_link(fit$link)
      tt <- terms(fit$formulas$formula.gamma)
      formula.2 <- formula(delete.response(tt))
      design_matrices <- model.matrix(formula.2, data = newdata)
      lps <- design_matrices %*% fit$coefs[[1]]
      pi <- link_fun(lps)
      pi
    }
  }else if(type == "probcure"){
    surv_fun <- get_surv(fit$dist)
    if(is.null(newdata)){
      pi <- link_fun(fit$coefs[[1]])
      lps <- t(matrix(unlist(fit$coefs[-1])))
      probtime <- pi / (pi + (1 - pi) * surv_fun(times, lp.k1 = lps[1], lp.k2 = lps[2], lp.k3= lps[3]))
      probcure <- data.frame(times = times, RS = probtime)
      probcure
    }else{
      tt <- terms(fit$formulas$formula.gamma)
      formula.2 <- formula(delete.response(tt))
      fit$formulas$formula.gamma <- formula.2
      design_matrices <- lapply(fit$formulas, function(x){
        if(length(x) > 0){
          return(model.matrix(x, data = newdata))
        }
      })

      lps <- vector("list", length(design_matrices))
      for(i in 1:length(lps)){
        if(length(design_matrices[[i]]) > 0){
          lps[[i]] <- design_matrices[[i]] %*% fit$coefs[[i]]
        }
      }

      #lps <- do.call(cbind, lps)
      probtime <- matrix(nrow = length(times), ncol = nrow(newdata))
      for(i in 1:nrow(newdata)){
        pi <- link_fun(lps[[1]][i,])
        probtime[, i] <- pi / (pi + (1 - pi) * surv_fun(times, lp.k1 = lps[[2]][i,], lp.k2 = lps[[3]][i,], lp.k3 = lps[[4]][i,]))
      }
      colnames(probtime) <- paste0("CP", 1:ncol(probtime))
      probcure <- cbind(times, probtime)
      probcure
    }
  }
}
