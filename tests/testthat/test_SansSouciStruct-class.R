test_that("Correctness of the constructor of SansSouciStruct", {
  s <- 100
  q <- 7
  m <- s * 2^q

  dirs <- c("bottom-up", "top-down")
  direction <- sample(dirs, 1)
  dir12 <- switch(direction,
    "bottom-up" = 1,
    "top-down" = 2
  )

  expect_warning(SansSouciDyadic(m,
    leaf_size = s, height = 8,
    direction = direction
  ))
  obj <- SansSouciDyadic(m, leaf_size = s, direction = direction)
  expect_s3_class(obj, "SansSouciStruct")

  dd <- dyadic_from_window_size(m, s, method = dir12)
  obj1 <- SansSouciStruct(dd$forest_structure, dd$leaf_list)
  expect_s3_class(obj, "SansSouciStruct")

  expect_identical(obj1, obj)
  rm(obj1)

  obj <- SansSouciDyadic(m, height = 8, direction = direction)
  expect_s3_class(obj, "SansSouciStruct")

  obj1 <- SansSouciStruct(dd$forest_structure, dd$leaf_list)
  expect_s3_class(obj, "SansSouciStruct")

  expect_identical(obj1, obj)
  rm(obj1)

  expect_equal(names(obj), c("input", "parameters", "output"))
  expect_null(obj$parameters)
  expect_null(obj$output)

  input <- obj$input
  expect_equal(names(input), c("struct", "leaves", "m"))
  expect_identical(input$struct, dd$forest_structure)
  expect_identical(input$leaves, dd$leaf_list)
  expect_length(input$leaves, 2^q)
  expect_length(input$struct, q + 1)

  expect_equal(n_hyp(obj), m)
  expect_equal(obj$input$m, m)
})


test_that("'fit.SansSouciStruct' reproduces the results of 'v_star'", {
  s <- 100
  q <- 7
  m <- s * 2^q

  obj <- SansSouciDyadic(m, leaf_size = s, direction = "top-down")

  K1 <- 8
  r <- 0.9
  m1 <- r * K1 * s
  barmu <- 3

  mu <- gen_mu_leaves(
    m = m, K1 = K1, d = r, grouped = TRUE,
    setting = "const",
    barmu = barmu, leaf_list = obj$input$leaves
  )
  mu_no <- gen_mu_noleaves(m = m, pi0 = 1, barmu = barmu)
  gb <- gauss_bloc(barmu, 10)

  obj$input$truth <- as.numeric(mu != 0)
  p_values <- gen_p_values(m = m, mu = mu, rho = 0)

  ord <- order(p_values)
  idxs <- round(c(
    seq(from = 1, to = 2 * m1, length = 30),
    seq(from = 2 * m1, to = m, length = 10)[-1]
  ))

  alpha <- 0.1

  # Oracle 1
  H0 <- which(mu == 0)
  V <- cumsum(ord %in% H0)
  V <- V[idxs]

  # Oracle 2
  res <- fit(obj, alpha = alpha, p_values = p_values, family = "Oracle")
  expect_s3_class(res, "SansSouci")
  v_mat <- predict(res, what = "FP", all = TRUE)
  VV <- v_mat$bound[idxs]
  expect_identical(VV, V)

  # Simes 1
  v_mat <- conf_curve_from_fam(
    p_values = p_values, ref_family = "Simes",
    param = alpha, what = "FP"
  )
  V <- v_mat$bound[idxs]

  # Simes 2
  res_simes <- fit(obj, alpha = alpha, p_values = p_values, family = "Simes")
  v_mat <- predict(res_simes, what = "FP", all = TRUE)
  VV <- v_mat$bound[idxs]
  expect_identical(VV, V)

  idxs1 <- head(idxs, 10)

  # DKWM 1 (tree)
  RES_DKWM <- fit(obj, alpha,
    p_values = p_values,
    family = "DKWM"
  )
  V <- sapply(idxs1, FUN = function(ii) {
    predict(RES_DKWM, S = seq_len(ii), what = "FP")
  })

  # DKWM 2 (tree)
  struct <- obj$input$struct
  leaves <- obj$input$leaves
  ZL <- structured_zetas(struct, leaves, structured_zeta_dkwm, p_values, alpha = alpha)
  VV <- sapply(idxs1, FUN = function(ii) {
    pred_pathway_structured(ord[1:ii],
                            forest_structure = struct,
                            zetas_list = ZL,
                            leaf_list = leaves
    )
  })
  expect_identical(VV, V)

  # DKWM 1 (part)
  RES_DKWM <- fit(obj, alpha,
    p_values = p_values,
    family = "DKWM", flavor = "partition"
  )
  V <- sapply(idxs1, FUN = function(ii) {
    predict(RES_DKWM, S = seq_len(ii), what = "FP")
  })

  # DKWM 2 (part)
  struct <- obj$input$struct
  CC <- struct[length(struct)]
  leaves <- obj$input$leaves
  ZL <- zetas_tree(CC, leaves, zeta.DKWM, p_values, alpha = alpha)

  VV <- sapply(idxs1, FUN = function(ii) {
    V.star(ord[1:ii],
      C = CC,
      ZL = ZL,
      leaf_list = leaves
    )
  })
  expect_identical(VV, V)
})
