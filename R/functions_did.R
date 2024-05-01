#' @export
generate_df_no_covariates <-
  function(
    N,
    T,
    T0
    ) {
      id <-
        matrix(
          rep(
            1:N,
            T
          ),
          N,
          T
        )

      time <-
        t(
          matrix(
            rep(
              1:T,
              N
            ),
            T,
            N
          )
        )

      a_i <-
        rep(stats::rnorm(
          n = N,
          mean = 5.4,
          sd = 0.15
        ),T)

      b_t <-
        rep(
          stats::rnorm(
            n = T,
            mean = 5.4 * 0.02,
            sd = 0.15
          ),
          N
        )

      e_it <- matrix(
        stats::rnorm(N * T,
          mean = 0,
          sd = 0.03
        ),
        N, T
      )

      U_i <-
        rep(stats::runif(
          n = N,
          min = 0,
          max =1
        ),T)

      y0_it <-
        matrix(
          a_i,
          N,
          T
        ) +
        t(
          matrix(
            b_t,
            T,
            N
          )
        ) +
        e_it

      g_i <-
        (
          U_i <=
          exp(-14 + 2.5 * a_i) /
          (1 + exp(-14 + 2.5 * a_i))
        )

      g_i <-
        matrix(
          g_i,
          N,
          T
        )

      d_it <-
        t(
          matrix(
            rep(
              c(
                rep(
                  0,
                  T0
                ),
                rep(
                  1,
                  T-T0
                )
              ),
              N
            ),
            T,
            N
          )
        ) * g_i

      tau_i <-
        rep(stats::rnorm(
          n = N,
          mean = 0.05,
          sd = 0.2
        ),T)

      tau_t <-
        t(
          matrix(
            rep(
              stats::runif(
                n = T,
                min = 0.9,
                max = 1.1
              ),
              N
            ),
            T,
            N
          )
        )

      tau_it <-
        matrix(
          tau_i,
          N,
          T
        ) *
        tau_t

      y_it <-
        y0_it +
        d_it *
        tau_it

      df <-
        tibble::tibble(
          id = as.vector(t(id)),
          time = as.vector(t(time)),
          a_i = as.vector(t(a_i)),
          g_i = as.vector(t(g_i)),
          d_it = as.vector(t(d_it)),
          y0_it = as.vector(t(y0_it)),
          y_it = as.vector(t(y_it)),
          U_i = as.vector(t(U_i)),
          tau_i = as.vector(t(tau_i)),
          tau_it = as.vector(t(tau_it))
        )
      return(df)
  }

#' @export
generate_df_multiperiod <-
  function(
    N, 
    T0, 
    T1, 
    T, 
    diff_trend = FALSE
  ) {
    id <- rep(1:N)
    a0 <- 5.4
    a_i <- 
      rnorm(
        N,
        mean = 0,
        sd = 0.15
      )
    U_i <- 
      runif(
        N, 
        0, 
        1
      )
    Ux_i <- 
      runif(
        N, 
        0, 
        1
      )
    x1_i <- (Ux_i >= 0.3)
    x2_i <- (Ux_i >= 0.7)
    tau_i <- 
      rnorm(
        N,
        mean = 0.1,
        sd = 0.2
      )

    group_i <- 0
    g1_i <- (
      U_i <= exp(-0.25) / (1 + exp(-0.25))
    )
    group_i <- 
      group_i + 
      g1_i * (T0 + 1)

    g2_i <- (
      (
        U_i > exp(-0.25) / (1 + exp(-0.25))
      ) & (
        U_i <= exp(-0.25 + 0.75 * (x1_i + 0)) / (1 + exp(-0.25 + 0.75 * (x1_i + 0)))
      )
    )
    group_i <- 
      group_i + 
      g2_i * (T0 + 2)

    g3_i <- (
      (
        U_i > exp(-0.25 + 0.75 * (x1_i + 0)) / (1 + exp(-0.25 + 0.75 * (x1_i + 0)))
      ) & (
        U_i <= exp(-0.25 + 0.75 * (x1_i + x2_i)) / (1 + exp(-0.25 + 0.75 * (x1_i + x2_i)))
      )
    )
    group_i <- 
      group_i + 
      g3_i * (T0 + 3)

    g4_i <- (
      (
        U_i > exp(-0.25 + 0.75 * (x1_i + x2_i)) / (1 + exp(-0.25 + 0.75 * (x1_i + x2_i)))
      ) & (
        U_i <= exp(-0.25 + 1 * (x1_i + x2_i)) / (1 + exp(-0.25 + 1 * (x1_i + x2_i)))
      )
    )
    group_i <- 
      group_i + 
      g4_i * (T0 + 4)

    for (
      t in seq(1, T)
    ) {
      time <- rep(t, N)
      if (diff_trend) {
        b_t <- (t / T) * (1 - x1_i - x2_i) - 0.2 * (t / T) * x1_i - 0.1 * (t / T) * x2_i
      } else {
        b_t <- t / T
      }
      e_it <- 
        rnorm(
            N,
            mean = 0,
            sd = 0.03
        )
      tau_t <- 
        runif(
          T,
          min = 0.9,
          max = 1.1
        ) 
      tau_it <- 
        tau_t * (
          abs(tau_i) * g1_i - 
          2.5 * abs(tau_i) * g2_i - 
          1.75 * abs(tau_i) * g3_i - 
          1 * abs(tau_i) * g4_i
        )

      d_it <- rep(0, N)
      if (t >= T0 + 1) {
        d1_it <- g1_i
      } else {
        d1_it <- rep(0, N)
      }
      if (t >= T0 + 2) {
        d2_it <- g2_i
      } else {
        d2_it <- rep(0, N)
      }
      if (t >= T0 + 3) {
        d3_it <- g3_i
      } else {
        d3_it <- rep(0, N)
      }
      if (t >= T0 + 4) {
        d4_it <- g4_i
      } else {
        d4_it <- rep(0, N)
      }

      d_it <- 
        d1_it + 
        d2_it + 
        d3_it + 
        d4_it

      y0_it <- 
        a0 + 
        (g1_i + g2_i + g3_i + g4_i) * (-a_i) + 
        (1 - g1_i - g2_i - g3_i - g4_i) * (a_i) + 
        b_t + 
        e_it
      y_it <- 
        y0_it + 
        d_it * tau_it

      df_period <- 
        tibble::tibble(
          id = id,
          x1_i = x1_i,
          x2_i = x2_i,
          group_i = group_i,
          g1_i = g1_i,
          g2_i = g2_i,
          g3_i = g3_i,
          g4_i = g4_i,
          time = time,
          tau_it = tau_it,
          d_it = d_it,
          d1_it = d1_it,
          d2_it = d2_it,
          d3_it = d3_it,
          d4_it = d4_it,
          y0_it = y0_it,
          y_it = y_it
        )
      if (
        t == 1
      ) {
        df_panel <- df_period
      } else {
        df_panel <- 
          rbind(
            df_panel, 
            df_period
          )
      }
    }
    return(df_panel)
  }

#' @export
generate_df_multiperiod_nyt <-
  function(
    N, 
    T0, 
    T1, 
    T, 
    diff_trend = FALSE
  ) {
    id <- rep(1:N)
    a0 <- 5.4
    a_i <- 
      rnorm(
        N,
        mean = 0,
        sd = 0.15
      )
    U_i <- 
      runif(
        N, 
        0, 
        1
      )
    Ux_i <- 
      runif(
        N, 
        0, 
        1
      )
    x1_i <- (Ux_i >= 0.3)
    x2_i <- (Ux_i >= 0.7)
    tau_i <- 
      rnorm(
        N,
        mean = 0.1,
        sd = 0.2
      )

    group_i <- 0
    g1_i <- (
      U_i <= exp(-0.25) / (1 + exp(-0.25))
    )
    group_i <- 
      group_i + 
      g1_i * (T0 + 1)

    g2_i <- (
      (
        U_i > exp(-0.25) / (1 + exp(-0.25))
      ) & (
        U_i <= exp(-0.25 + 0.75 * (x1_i + 0)) / (1 + exp(-0.25 + 0.75 * (x1_i + 0)))
      )
    )
    group_i <- 
      group_i + 
      g2_i * (T0 + 2)

    g3_i <- (
      (
        U_i > exp(-0.25 + 0.75 * (x1_i + 0)) / (1 + exp(-0.25 + 0.75 * (x1_i + 0)))
      ) & (
        U_i <= exp(-0.25 + 0.75 * (x1_i + x2_i)) / (1 + exp(-0.25 + 0.75 * (x1_i + x2_i)))
      )
    )
    group_i <- 
      group_i + 
      g3_i * (T0 + 3)

    g4_i <- (
      (
        U_i > exp(-0.25 + 0.75 * (x1_i + x2_i)) / (1 + exp(-0.25 + 0.75 * (x1_i + x2_i)))
      ) & (
        U_i <= exp(-0.25 + 1 * (x1_i + x2_i)) / (1 + exp(-0.25 + 1 * (x1_i + x2_i)))
      )
    )
    group_i <- 
      group_i + 
      g4_i * (T0 + 4)

    group_i <- 
      group_i + 
      (1 - g1_i - g2_i - g3_i - g4_i) * T

    for (
      t in seq(1, T)
    ) {
      time <- rep(t, N)
      if (diff_trend) {
        b_t <- 
          (t / T) * (1 - x1_i - x2_i) - 
          0.2 * (t / T) * x1_i - 
          0.1 * (t / T) * x2_i
      } else {
        b_t <- t / T
      }
      e_it <- 
        rnorm(
          N,
          mean = 0,
          sd = 0.03
        )
      tau_t <- 
        runif(
          T,
          min = 0.9,
          max = 1.1
        )
      tau_it <- 
        tau_t * (
          abs(tau_i) * g1_i - 
          2.5 * abs(tau_i) * g2_i - 
          1.75 * abs(tau_i) * g3_i - 
          1 * abs(tau_i) * g4_i
        )

      d_it <- rep(0, N)
      if (t >= T0 + 1) {
        d1_it <- g1_i
      } else {
        d1_it <- rep(0, N)
      }
      if (t >= T0 + 2) {
        d2_it <- g2_i
      } else {
        d2_it <- rep(0, N)
      }
      if (t >= T0 + 3) {
        d3_it <- g3_i
      } else {
        d3_it <- rep(0, N)
      }
      if (t >= T0 + 4) {
        d4_it <- g4_i
      } else {
        d4_it <- rep(0, N)
      }

      d_it <- 
        d1_it + 
        d2_it + 
        d3_it + 
        d4_it
      if (t == T) {
        d_it <- 1
      }

      y0_it <- 
        a0 + 
        (g1_i + g2_i + g3_i + g4_i) * (-a_i) + 
        (1 - g1_i - g2_i - g3_i - g4_i) * (a_i) + 
        b_t + 
        e_it
      y_it <- 
        y0_it + 
        d_it * tau_it

      df_period <- 
        tibble::tibble(
          id = id,
          x1_i = x1_i,
          x2_i = x2_i,
          group_i = group_i,
          g1_i = g1_i,
          g2_i = g2_i,
          g3_i = g3_i,
          g4_i = g4_i,
          time = time,
          tau_it = tau_it,
          d_it = d_it,
          d1_it = d1_it,
          d2_it = d2_it,
          d3_it = d3_it,
          d4_it = d4_it,
          y0_it = y0_it,
          y_it = y_it
        )
      if (t == 1) {
        df_panel <- df_period
      } else {
        df_panel <- 
          rbind(
            df_panel, 
            df_period
          )
      }
    }
    return(df_panel)
  }
