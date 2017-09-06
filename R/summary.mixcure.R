
#' Summary of parametric cure models
#'
#' This function generates a summary of the fitted cure model.
#'
#' @param formula An object of type \code{CureModel}
#' @return An object of class \code{summary.CureModel}.
#' @export

summary.CureModel <- function(fit){
  se <- sqrt(diag(fit$cov))
  tval <- unlist(fit$coefs) / se
  coefs <- unlist(fit$coefs)
  results <- lapply(c("gamma", "k1", "k2", "k3"), function(x){
    these <- grepl(x, names(se))
    coef <- coefs[these]
    TAB <- cbind(Estimate = coef,
                 StdErr = se[these],
                 t.value = tval[these],
                 p.value = 2*pt(-abs(tval[these]), df = fit$df))

    res <- list(call = fit$formula[paste0("formula.", x)],
                coefficients = TAB)
  })
  results <- list(results)
  results$type <- fit$type
  results$link <- fit$link
  results$dist <- fit$dist
  results$ML <- fit$ML
  class(results) <- "summary.CureModel"
  results
}

print.summary.CureModel <- function(x)
{
  for(i in 1:length(x[[1]])){
    print_call <- paste0("Call - ", names(x[[1]][[i]]$call),":")
    cat(print_call, "\n")
    print(x[[1]][[i]]$call[[1]])
    #    cat("\n")
    if(i == length(x[[1]])){
      printCoefmat(x[[1]][[i]]$coefficients, P.value = TRUE, has.Pvalue = TRUE)
    }else{
      printCoefmat(x[[1]][[i]]$coefficients)
    }
    cat("\n")
  }
  cat("Type =", x$type, "\n")
  cat("Distribution =", x$dist, "\n")
  cat("Link =", x$link, "\n")
  cat("LogLik(model) =", x$ML, "\n")

}
