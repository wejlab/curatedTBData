returned_resources <- curatedTBData(c("GSE107104", "GSE19435", "GSE19443"),
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
            expect_error(boxplotTBSig(re1_noname, gset = "Sweeney_OD_3",
                                      annotationColName = "TBStatus"))
            expect_error(boxplotTBSig(re1_empty_name, gset = "Sweeney_OD_3",
                                      annotationColName = "TBStatus"))
          })
