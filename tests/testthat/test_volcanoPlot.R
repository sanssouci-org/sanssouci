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
  expect_equal(attr(vp, "selection")$sel12, seq_len(m))     ## no active filter: all genes selected

  vp <- volcanoPlot(cal, p = 1, q = 1, r = 0, show_signed_selections = TRUE)
  expect_is(vp, "ggplot")
  expect_equal(attr(vp, "selection")$sel12, seq_len(m))     ## no active filter: all genes selected

  vp <- volcanoPlot(cal, p = 1, q = 1, r = Inf)
  expect_equal(length(vp$selection), 0)      ## too stringent filtering on fold change

  vp <- volcanoPlot(cal, p = 0, q = 1, r = 0)
  expect_equal(length(vp$selection), 0)      ## too stringent filtering on p-value

  vp <- volcanoPlot(cal, p = 1, q = 0, r = 0)
  expect_equal(length(vp$selection), 0)     ## too stringent filtering on q-value

  vp <- volcanoPlot(cal, p = 1, q = 1, r = 0)
  expect_null(vp$selection)
  expect_is(vp, "ggplot")

  vp_build <- ggplot2::ggplot_build(vp)

  # Test 1 : vérifier le format de text_hover
  point_data <- vp_build$data[[5]]

  expect_match(
    point_data$text,
    "^feature: \\d+<br>logFC: -?\\d+(\\.\\d+)?([eE][+-]?\\d+)?<br>p-value: \\d+(\\.\\d+)?([eE][+-]?\\d+)?$"
  )

  # Test 2 : check mapping text and colone in dataframe
  expect_identical(
    rlang::get_expr(vp$layers[[5]]$mapping$text),
    quote(.data$text_hover)
  )


  ## filtering both on p-values and q-values
  expect_warning(vp <- volcanoPlot(cal, p = 0.01, q = 0.05, r = 0,
                                   ylim = c(0, 6)))

})

test_that("Convert volcanoplot with ggploty", {
  m <- 500
  pi0 <- 0.5
  m1 <- m-m*pi0
  SNR <- 5*(runif(m1)-0.5)
  obj <- SansSouciSim(m = m, rho = 0.4, n = 100,
                      pi0 = pi0, SNR = SNR, prob = 0.5)
  alpha <- 0.2
  cal <- fit(obj, alpha = alpha, B = 1e2, family="Simes")

  vp <- volcanoPlot(cal, p = 1, q = 1, r = 0)

  vpi <- plotly::ggplotly(vp)
  expect_s3_class(vpi, "plotly")

  ## expect reproduce ggplot2 plot with all information
  feature_label = "feature"
  sel <- attr(vp, "selection")
  vp <- plotly::ggplotly(vp +
                           ggplot2::labs(
                             y = "p-value (-log<sub>10</sub> scale)"
                           ),
                         tooltip = "text"
  ) |>
    plotly::layout(
      margin = list(t = 60, b = 70)
    )

  txt_right <- sprintf(
    "%d %s\nTP \u2265 %d ; FDP \u2264 %.2f",
    sel$n1, pluralize(word = feature_label, n = sel$n1), sel$TP1, sel$FDP1
  )

  txt_left <- sprintf(
    "%d %s\nTP \u2265 %d ; FDP \u2264 %.2f",
    sel$n2, pluralize(word = feature_label, n = sel$n2), sel$TP2, sel$FDP2
  )
  vp <- vp |>
    plotly::layout(
      annotations = list(
        list(
          x = 1,
          y = 1,
          xref = "paper",
          yref = "paper",
          text = txt_right,
          showarrow = FALSE,
          xanchor = "right",
          yanchor = "top",
          align = "right"
        ),
        list(
          x = 0,
          y = 1,
          xref = "paper",
          yref = "paper",
          text = txt_left,
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "top",
          align = "left"
        )
      )
    )
  expect_s3_class(vp, "plotly")

})


