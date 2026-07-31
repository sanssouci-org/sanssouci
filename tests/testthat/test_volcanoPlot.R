context("Volcano plots")

test_that("Vanilla test for 'volcanoPlot'", {
  m <- 500
  pi0 <- 0.5
  m1 <- m-m*pi0
  SNR <- 5*(runif(m1)-0.5)
  obj <- SansSouciSim(m = m, rho = 0.4, n = 100,
                      pi0 = pi0, SNR = SNR, prob = 0.5)
  alpha <- 0.2
  cal <- fit(obj, alpha = alpha, B = 1e2, family="Simes")

  vp <- volcanoPlot(cal, p = 1, q = 1, r = 0)
  expect_is(vp, "ggplot")
  expect_equal(attr(vp, "selection"), seq_len(m))     ## no active filter: all genes selected

  vp <- volcanoPlot(cal, p = 1, q = 1, r = 0, show_signed_selections = TRUE)
  expect_is(vp, "ggplot")
  expect_equal(attr(vp, "selection"), seq_len(m))     ## no active filter: all genes selected

  vp <- volcanoPlot(cal, p = 1, q = 1, r = Inf)
  expect_equal(length(vp$selection), 0)      ## too stringent filtering on fold change

  vp <- volcanoPlot(cal, p = 0, q = 1, r = 0)
  expect_equal(length(vp$selection), 0)      ## too stringent filtering on p-value

  vp <- volcanoPlot(cal, p = 1, q = 0, r = 0)
  expect_equal(length(vp$selection), 0)     ## too stringent filtering on q-value

  vp <- volcanoPlot(cal, p = 1, q = 1, r = 0)
  expect_null(vp$selection)
  expect_is(vp, "ggplot")


  expect_match(
    vp@data$text_hover[1],
    "^feature: \\d+<br>logFC: -?\\d+(\\.\\d+)?([eE][+-]?\\d+)?<br>p-value: \\d+(\\.\\d+)?([eE][+-]?\\d+)?$"
  )
  expect_identical(
    rlang::get_expr(vp@layers$geom_point$mapping$text),
    quote(.data$text_hover)
  )

  ## filtering both on p-values and q-values
  expect_warning(vp <- volcanoPlot(cal, p = 0.01, q = 0.05, r = 0,
                                   ylim = c(0, 6)))

  vp <- volcanoPlot(cal, p = 0.5, q = 1, r = 0, interactive = TRUE)
  expect_null(vp$selection)
  expect_is(vp, "plotly")
})
