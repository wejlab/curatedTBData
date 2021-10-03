returned_resources <- curatedTBData(c("GSE107104", "GSE19435"),
                                    dry.run = FALSE, curated.only = TRUE) |>
    suppressWarnings()
mysignatures <- list(Sweeney_OD_3 = c("DUSP3", "GBP5", "KLF2"),
                     Suliman_4 = c("ANKRD22", "C1QC", "OSBPL10", "TRAV27"))
sig_names <- names(mysignatures)
re1 <- lapply(returned_resources, function(x)
    subset_curatedTBData(x, annotationColName = "TBStatus",
                         annotationCondition = c("Control", "PTB")))
re2 <- lapply(re1, function(x)
    TBSignatureProfiler::runTBsigProfiler(input = x,
                                          useAssay = 1,
                                          signatures = mysignatures,
                                          algorithm = "ssGSEA",
                                          update_genes = FALSE))

test_that("Input list does not have unique name for each element.
          input list has empty name for one or more elements", {
              re1_noname <- re1
              names(re1_noname) <- NULL
              re1_empty_name <- re1
              names(re1_empty_name)[1:2] <- ""
              expect_error(combine_auc(re1_noname,
                                       annotationColName = "TBStatus",
                                       signatureColNames = sig_names))
              expect_error(combine_auc(re1_empty_name,
                                       annotationColName = "TBStatus",
                                       signatureColNames = sig_names))
})

test_that("Return error when \"annotationColName\" is not found within study", {
    expect_error(combine_auc(re2, annotationColName = "fakeAnnotation",
                             signatureColNames = sig_names)) |>
        suppressWarnings()
})

test_that("The number of levels for annotation column should be exactly 2", {
    re3 <- lapply(re2, function(x) {
        n <- dim(x)[2]
        x$fakeAnnotation1 <- letters[1]
        x$fakeAnnotation2 <- rep(letters[1:3], ceiling(n / 3))[1:n]
        x
    })
    expect_error(combine_auc(re3, annotationColName = "fakeAnnotation1",
                             signatureColNames = sig_names) |>
                     suppressWarnings())
    expect_error(combine_auc(re3, annotationColName = "fakeAnnotation2",
                             signatureColNames = sig_names) |>
                     suppressWarnings())
})

test_that("Return error when all \"signatureColNames\" are not found
          within each study", {
    expect_error(combine_auc(re2, annotationColName = "TBStatus",
                            signatureColNames = "fakeName"))
})

test_that("Message when some \"signatureColNames\" are not found
          within study", {
    expect_message(combine_auc(re2, annotationColName = "TBStatus",
                               signatureColNames = c(sig_names, "fakeName")))
})

test_that("Message for not computing Bootstrap Confidence Interval", {
    expect_message(combine_auc(re2, annotationColName = "TBStatus",
                               signatureColNames = sig_names))
})

test_that("Message for when constant score is found for signatures", {
    re3 <- lapply(re2, function(x) {
        x$fakeName <- 0.5
        x
    })
    expect_message(combine_auc(re3, annotationColName = "TBStatus",
                               signatureColNames = c(sig_names, "fakeName")))
})

test_that("Return class is data.frame", {
    out <- combine_auc(re2, annotationColName = "TBStatus",
                       signatureColNames = sig_names)
    expect_s3_class(out, "data.frame")
})
