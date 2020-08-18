context("Tree functions")

test_that("Vanilla test for tree functions", {
    m <- 6
    check_tree <- function(dd) {
        expect_true(class(dd) == "list")
        expect_length(dd, 2)
        expect_named(dd)
        expect_identical(names(dd), c("leaf_list", "C"))
    }

    dd <- dyadic.from.window.size(m, s = 2, method = 2)
    check_tree(dd)
    
    dd <- dyadic.from.height(m, H = 3, method = 2)
    check_tree(dd)
    
    dd <- dyadic.from.max.height(m, method=2)
    check_tree(dd)
})