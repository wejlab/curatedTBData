returned_resources <- curatedTBData(c("GSE107104", "GSE19435"),
                                    dry.run = FALSE, curated.only = TRUE) |>
    suppressWarnings()
mysignatures <- list(Sweeney_OD_3 = c("DUSP3", "GBP5", "KLF2"))
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
            expect_error(boxplotTBSig(re1_noname,
                                      annotationColName = "TBStatus",
                                      signatureColNames = "Sweeney_OD_3"))
            expect_error(boxplotTBSig(re1_empty_name,
                                      annotationColName = "TBStatus",
                                      signatureColNames = "Sweeney_OD_3"))
})

test_that("Signature name is not found from the entire input", {
    expect_error(boxplotTBSig(re2, annotationColName = "TBStatus",
                              signatureColNames = "fakeGeneSignatures"))

})

test_that("Annotation name is not found from the entire input", {
    expect_error(boxplotTBSig(re2, annotationColName = "fakeAnnotationName",
                              signatureColNames = "Sweeney_OD_3"))

})

test_that("The number of condition levels are greater than 9", {
    fake_level <- letters[1:10]
    re3 <- lapply(re2, function(x) {
        n <- dim(x)[2]
        x$myAnnotation <- c(fake_level, sample(fake_level, n - 10, TRUE))
        x
    })
    expect_message(boxplotTBSig(re3, annotationColName = "myAnnotation",
                                signatureColNames = "Sweeney_OD_3"))
})

test_that("Return type is a list of gtable", {
    p <- boxplotTBSig(re2, annotationColName = "TBStatus",
                      signatureColNames = "Sweeney_OD_3")
    expect_type(p, "list")
    expect_s3_class(p, "gtable")
})
