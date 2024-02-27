#'@export
return_dgp_cluster <-
  function(
    N,
    N_C
  ) {
    X <-
      rnorm(
        N * N_C
      )
    C <-
      rep(
        1:N_C,N
      )
    e_c <-
      rnorm(
        N_C
      )
    lmd <-
      rnorm(
        N * N_C,
        mean = 2
      )^2
    tau <-
      seq(
        -2,
        2,
        length = N_C
      )
    e <-
      lmd * e_c[C] +
      rnorm(
        N * N_C,
        sd = 0.1
      )
    Y <-
      0.5 +
      X * tau[C] +
      e

    return(
      data.frame(
        Y = Y,
        X = X,
        C = C
      )
    )
}
