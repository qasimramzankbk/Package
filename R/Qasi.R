
#' Fahrenheit conversion
#'
#' Convert degrees Fahrenheit temperatures to degrees Celsius
#' @param F_temp The temperature in degrees Fahrenheit
#' @return The temperature in degrees Celsius
#' @examples
#' temp1 <- F_to_C(50);
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export
F_to_C <- function(F_temp){
  C_temp <- (F_temp - 32) * 5/9;
  return(C_temp);
}



#' Celsius conversion
#'
#' Convert degrees Celsius temperatures to degrees Fahrenheit
#' @param C_temp The temperature in degrees Celsius
#' @return The temperature in degrees Fahrenheit
#' @examples
#' temp1 <- C_to_F(22);
#' temp2 <- C_to_F( c(-2, 12, 23) );
#' @export
C_to_F <- function(C_temp){
  F_temp <- (C_temp * 9/5) + 32;
  return(F_temp);
}




#' Quantile Function for Modified Weibull Distribution (MWD)
#'
#' This function generates random observations by calculating the Quantile Function for the Modified Weibull Distribution (MWD), given the shape, scale, and rate parameters.
#' The function uses the inverse transformation method to generate random samples from the specified distribution.
#'
#' @param n Integer. The number of random observations to generate.
#' @param alpha Numeric. The shape parameter of the Modified Weibull Distribution (MWD). Must be greater than 1.
#' @param beta Numeric. The scale parameter of the Modified Weibull Distribution (MWD). Must be positive.
#' @param lambda Numeric. The rate parameter of the Modified Weibull Distribution (MWD). Must be positive.
#'
#' @details
#' The quantile function is the inverse of the cumulative distribution function (CDF) and is often used in simulation studies and reliability analysis.
#' This function generates random numbers by calculating the quantile function of the Modified Weibull Distribution (MWD), which is particularly useful for data modeling with a variety of real-world applications such as survival analysis, reliability engineering, and risk analysis.
#' The method used to calculate quantiles is based on the formula for the quantile function:
#' \deqn{Q(p) = \left( \frac{-log(1 - \frac{log(p(\alpha - 1) + 1)}{log(\alpha)})}{\lambda} \right)^{1/\beta}}
#' where \eqn{p} is a probability value.
#'
#' @return
#' A numeric vector of length \code{n} containing generated observations from the Modified Weibull Distribution (MWD).
#'
#' @examples
#' # Generate 4 random samples from MWD with specified parameters
#' temp1 <- rqasi(4, 0.2, 0.3, 0.4)
#'
#' # Generate 40 random samples from MWD with different parameters
#' temp2 <- rqasi(40, 2, 0.3, 1.4)
#'
#' @seealso \code{\link{rweibull}}, \code{\link{runif}}
#'
#' @references
#' Ramzan, Q., Amin, M., & Faisal, M. (2022). Bayesian inference for modified Weibull distribution under simple step‐stress model based on type‐I censoring. \emph{Quality and Reliability Engineering International}, 38(2), 757-779.
#'
#' @note
#' Ensure that the parameters \code{alpha}, \code{beta}, and \code{lambda} are chosen appropriately to avoid numerical issues.
#' @export

rqasi <- function(n, alpha, beta, lambda){
  quantile_function <- function(p, alpha, beta, lambda) {
    Q <- function(p) {
      ((-log(1 - (log(p * (alpha - 1) + 1) / log(alpha))) / lambda)^(1 / beta))
    }
    return(Q(p))
  }
  vectorized_q <- Vectorize(quantile_function, vectorize.args = c("p"))
  p <- runif(n)  # Generate n uniform random numbers between 0 and 1
  samples <- vectorized_q(p, alpha, beta, lambda)
  return(samples)
}





#' Probability Density Function for Qasi Distribution
#'
#' This function calculates the Probability Density Function (PDF) for the Qasi distribution, given the shape, scale, and rate parameters.
#' The Qasi distribution is a flexible probability model used in reliability analysis and survival studies, suitable for modeling real-world data involving failure times and lifetime distributions.
#'
#' @param t Numeric. The time or variable for which the PDF is to be computed.
#' @param alpha Numeric. The shape parameter of the Qasi Distribution. Must be greater than 1.
#' @param beta Numeric. The scale parameter of the Qasi Distribution. Must be positive.
#' @param lambda Numeric. The rate parameter of the Qasi Distribution. Must be positive.
#'
#' @details
#' The probability density function describes the likelihood of a random variable taking a specific value. In the Qasi distribution, the PDF is useful for estimating failure rates and assessing reliability in various applications.
#' The formula for the Qasi distribution PDF is given by:
#' \deqn{f(t) = 4\alpha\theta^{2} t e^{-2\theta t }\bigg[1-\left(1+2\theta t \right)e^{-2\theta t }\bigg]^{\alpha-1}}
#' where \eqn{t} is the time or variable of interest.
#'
#' @return
#' A numeric value representing the probability density for the given parameters at the specified \code{t}.
#'
#' @examples
#' # Calculate the PDF for a given time point with specified parameters
#' pdf_value <- pdf_qasi(2, 1.5, 0.8, 0.3)
#'
#' @seealso \code{\link{dweibull}}, \code{\link{dexp}}
#'
#' @references
#' Ramzan, Q., Amin, M., & Faisal, M. (2022). Bayesian inference for modified Weibull distribution under simple step‐stress model based on type‐I censoring. \emph{Quality and Reliability Engineering International}, 38(2), 757-779.
#'
#' @note
#' Ensure that the parameters \code{alpha}, \code{beta}, and \code{lambda} are chosen appropriately to avoid numerical issues or undefined behavior.
#' @export

dqasi <- function(t, beta, zeta, phi) {
  # Calculate the PDF based on the given formula
  pdf <- beta * t^(zeta - 1) * (zeta + phi * t) *
    exp(-beta * t^zeta * exp(phi * t) + phi * t)

  return(pdf)
}

