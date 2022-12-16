
Jacobian_by_diff <- function(f, x, xdata, eps)
{
  k = length(x)
  Y = f(x, xdata)
  dYdx = Reduce(function(p, i)
  {
    e = numeric(length(x))
    e[i] = max(1, abs(x[i])) * eps
    Yi = f(x + e, xdata)
    dYi = (Yi - Y) / e[i]
    return(cbind(p, dYi))
  },
  1:k, init = matrix(nrow = nrow(xdata), ncol = 0))
  f_name = deparse(substitute(f))
  x_names = names(x)
  if (is.null(x_names))
    x_names = 1:length(x)

  colnames(dYdx) <- paste0(f_name, '_', x_names)
  return (dYdx)
}

gradient_by_diff <- function(f, x, xdata, ydata=NULL,eps=1e-5)
{
  if (is.null(ydata))
    ydata=0
  k = length(x)
  Y = f(x, xdata)
  dYdx = Reduce(function(p, i)
  {
    e = numeric(length(x))
    e[i] = max(1, abs(x[i])) * eps
    Yi = f(x + e, xdata)
    dYi = (Yi - Y) / e[i]
    return(cbind(p, dYi))
  },
  1:k, init = matrix(nrow = nrow(xdata), ncol = 0))
  d_f__d_x= Conj(t(dYdx)) %*% (Y - ydata)


  f_name = deparse(substitute(f))
  x_names = names(x)
  if (is.null(x_names))
    x_names = 1:length(x)



  rownames(d_f__d_x) <- paste0(f_name, '_', x_names)
  return (d_f__d_x)
}


get_stochastic_indexes <-
  function(number_of_samples, number_of_chunks)
  {
    random_indexes = sample.int(number_of_samples)
    return(lapply(1:number_of_chunks, function (i)
    {
      i_start = 1 + (i - 1) / number_of_chunks * number_of_samples
      i_end = i / number_of_chunks * number_of_samples
      return(random_indexes[i_start:i_end])
    }))
  }







#' Stochastic Descent with ADAM step learning
#' Minimizes the square sum of f(x,xdata), taking samples of xdata for calculating the steps
#'
#' @param f  function to be optimized
#' @param x0 initial guess of the multidimentional parameter to be optimized
#' @param xdata  independent possibly multivariate variable where f is applied
#' @param ydata  rowise monovariate variable given for each xdata
#' @param Jacobian optional, function that produces the Jacobian of f
#' @param number_of_chunks number of samples from the data
#' @param maxiter maximum number of iteration
#' @param return_Ft  whether the function values are returned
#' @param tol_gradient termination tolerance for the inf norm of the gradient
#' @param tol_x  termination tolerance for changes in x
#' @param alpha  gradient descent step size
#' @param beta1  exponential decay rates for the mean estimates
#' @param beta2  exponential decay rates for the variance estimates
#' @param eps  minimum std estimate
#' @param eps_J epsilon for Jacobian
#'
#' @return a list containing the optimal x and other things
#' @export
#'
#' @examples
#'
#' # It applies the algorithm described in https://arxiv.org/pdf/1412.6980.pdf
#'
#' sum_of_Normals <- function(x, X){
#'   return(x["A1"] * exp(-(X[,"x1"] - x["mean1"])^2/x["sd1"]^2) +
#'          x["A2"] * exp(-(X[,"x2"] - x["mean2"])^2/x["sd2"]^2))
#' }
#'
#' x = c(A1 = 12.3,mean1 = 24, sd1 = 17.3, A2 = 0.23,mean2 = 21.2 ,sd2 = 5)
#' X=cbind(runif(1e5) * 100 - 20, runif(1e5) * 100 - 20)
#' colnames(X)<-c('x1','x2')
#'
#' Yn = sum_of_Normals(x,X)+rnorm(1e5)*0.2
#'
#'
#' opt=nlsqrsgd(f =sum_of_Normals, x0 = x, xdata = X, ydata = Yn,number_of_chunks = 100)
#'
nlsqrsgd <- function(f,
                     x0,
                     xdata,
                     ydata,
                     number_of_chunks,
                     maxiter = 10000,
                     is_complex = FALSE,
                     tol_gradient = 1e-5,
                     tol_x = 1e-5,
                     return_Ft = FALSE,
                     gradient = NULL,
                     Jacobian = NULL,
                     use_gradient = FALSE,
                     alpha = 0.001,
                     beta1 = 0.9,
                     beta2 = 0.999,
                     eps = 1e-8,
                     eps_J = 1e-5)
{
  if (is.complex(x0))
    is_complex = TRUE
  if (!is.matrix(ydata))
    ydata = matrix(ydata, ncol = 1)
  stopifnot("same number of rows" = nrow(xdata) == nrow(ydata))

  if (!is.null(gradient))
    use_gradient = TRUE

  if ((use_gradient) && (is.null(gradient)))
    gradient <- function(x,xdata,ydata)
      gradient_by_diff(f = f,
                       x = x,
                       xdata = xdata,
                       ydata = ydata,
                       eps = eps_J)

  if ((!use_gradient) && (is.null(Jacobian)))

    Jacobian <-
    function(x, xdata)
      Jacobian_by_diff(
        f = f,
        x = x,
        xdata = xdata,
        eps = eps_J
      )


    k = length(x0)


  initial_step <- function(x0)
  {
    if (is_complex)
    {
    mt = complex(k)
    vt =  complex(k)
    ca = list(
      xt = x0,
      xp = x0,
      mt = mt,
      vt = vt,
      sqrsum_tot = 0,
      gradient_tot = complex(k),
      done = FALSE
    )
    return(ca)
    }
    else
    {
      mt = numeric(k)
      vt =  numeric(k)
      ca = list(
        xt = x0,
        xp = x0,
        mt = mt,
        vt = vt,
        sqrsum_tot = 0,
        gradient_tot = numeric(k),
        done = FALSE)
      return(ca)

    }

  }



  n = nrow(xdata)
  stc_ind = get_stochastic_indexes(number_of_samples = n, number_of_chunks = number_of_chunks)
  Reduce(
    f = function(pr, t) {
      if ((pr$done))
        return(pr)
      ca = list(
        xp = pr$xp,
        nsteps = t,
        done = FALSE
      )
      ii = stc_ind[[(t - 1) %% number_of_chunks + 1]]
      Ft = f(pr$xt, xdata[ii,])
      yt = ydata[ii,]
      ca$sqrsum_t = Re(sum(Conj(Ft - yt)*(Ft - yt)))
      if (use_gradient)
        ca$gradient_t = gradient(pr$xt, xdata[ii,],yt)
      else
      {
        Jt = Jacobian(pr$xt, xdata[ii,])
        ca$gradient_t = Conj(t(Jt)) %*% (Ft - yt)
      }
      ca$sqrsum_tot = pr$sqrsum_tot + ca$sqrsum_t
      ca$gradient_tot = pr$gradient_tot + ca$gradient_t
      ca$mt = beta1 * pr$mt + (1 - beta1) * ca$gradient_t
      ca$vt = beta2 * pr$vt + (1 - beta2) * ca$gradient_t ^ 2
      alphat = alpha * sqrt(1 - beta2 ^ t) / (1 - beta1 ^ t)
      ca$xt = pr$xt
      ca$xt[seq_along(ca$xt)] = ca$xt- alphat * ca$mt / (sqrt(ca$vt) + eps)
      names(ca$xt) <- names(pr$xt)


      if (return_Ft)
      {
        ca$Ft_tot = pr$Ft_tot
        ca$Ft_tot[ii] = Ft
      }
      if (t %% number_of_chunks != 0)
        return(ca)
      else if (max(abs(Re(ca$gradient_tot))) < tol_gradient)
      {
        ca$stopcondition = list(reason = 'small gradient norm', gradient_norm =
                                  max(abs(ca$gradient_tot)))
        ca$done = TRUE
        return(ca)
      }
      else if (max(abs(ca$xp - ca$xt)) < tol_x)
      {
        ca$stopcondition = list(reason = 'no progress',
                                max_delta_x = max(abs(ca$xp - ca$xt)))
        ca$done = TRUE
        return(ca)
      }
      else
      {
        ca$xp = ca$xt
        ca$sqrsum_tot = 0
        if (t < maxiter)
          ca$gradient_tot = ca$gradient_tot - ca$gradient_tot
        return(ca)
      }
    },
    x = 1:maxiter,
    init = initial_step(x0 = x0)
  )

}
