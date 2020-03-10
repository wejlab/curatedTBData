#' calculate molecular degree of perturbation scores for selected sample

setGeneric("MDP", function(theObject, gset, ...) standardGeneric("MDP"))

setMethod("MDP", signature (theObject = "SummarizedExperiment", gset = "NULL"),
          function(theObject, gset = NULL, assay_name = 1){
   # theObject data after normalization
   counts <- assays(theObject)[[assay_name]]
   # select the reference group, calculate mean and sd
   theObject_ref <- assays(theObject[,theObject$TBStatus == "Control"])[[assay_name]]
   mean_ref <- apply(theObject_ref, 1, mean)
   sd_ref <- apply(theObject_ref, 1, sd)

   # calculate gMDP for each gene i and sample s
   gMDP <- apply(counts, 2, function(x){(x-mean_ref)/sd_ref})

   # replace small gMDP values, set abs(gMDP) < 2 equals to 0
   gMDP[abs(gMDP)<2] <- 0

   # Include the gMDP assays in the existing SummarizedExperiment Object
   assays(theObject)[["gMDP"]] <- gMDP

   # show available classes in datasets
   TBStatus <- unique(colData(theObject)$TBStatus %>% as.character())

   # calculate the average gMDP score wihtin each class and sort from max to minimal
   gMDP_avg <- lapply(1:length(TBStatus), function(x){
     gMDP_class <- assays(theObject[,theObject$TBStatus == TBStatus[x]])[["gMDP"]]
     gMDP_class_avg <- apply(gMDP_class, 1, mean) %>% sort(decreasing = T)
     gMDP_class_avg
   })
   names(gMDP_avg) <- TBStatus

   # Use the original method in the paper
   # Identify genes with top 25% gMDP_avg within each class
   num_25 <- floor(nrow(gMDP)*0.25)
   gMDP_avg_gene <- lapply(gMDP_avg, function(x) x[1:num_25] %>% names())

   # claculate MDP of individual sample (sMDP)
   sMDP <- lapply(1:length(TBStatus), function(x,gMDP_avg_gene){
     gMDP_class <- assays(theObject[,theObject$TBStatus == TBStatus[x]])[["gMDP"]]
     sMDP_class <- gMDP_class[which(row.names(gMDP_class) %in% gMDP_avg_gene[[x]]),]
     sMDP_sample <- apply(sMDP_class, 2, mean)
     sMDP_sample
   }, gMDP_avg_gene)
   names(sMDP) <- TBStatus

   sMDP_data <- do.call(rbind, lapply(1:length(sMDP), function(x) {
     data.frame(Sample=names(sMDP[[x]]), sMDP = as.vector(sMDP[[x]]), TBStatus = names(sMDP[x]))
   }))

   return(sMDP_data)

})

setMethod("MDP", signature (theObject = "MultiAssayExperiment", gset = "NULL"),
          function(theObject, gset = NULL, experiment_type = "assay_reduce", assay_name = 1){

            # Create SummarizedExperiment
            col_data = colData(theObject)
            if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
              index <- sapply(1:length(colnames(theObject[[experiment_type]])), function (i)
                which(row.names(col_data) %in% colnames(theObject[[experiment_type]])[i]))

              col_data <- col_data[index,]
            }

            theObject <- SummarizedExperiment::SummarizedExperiment(assays = list(counts= as.matrix(theObject[[experiment_type]])), colData = col_data)

            # Everything follows the summarized Experiment
            counts <- assays(theObject)[[assay_name]]
            # select the reference group, calculate mean and sd
            theObject_ref <- assays(theObject[,theObject$TBStatus == "Control"])[[assay_name]]
            mean_ref <- apply(theObject_ref, 1, mean)
            sd_ref <- apply(theObject_ref, 1, sd)

            # calculate gMDP for each gene i and sample s
            gMDP <- apply(counts, 2, function(x){(x-mean_ref)/sd_ref})

            # replace small gMDP values, set abs(gMDP) < 2 equals to 0
            gMDP[abs(gMDP)<2] <- 0

            # Include the gMDP assays in the existing SummarizedExperiment Object
            assays(theObject)[["gMDP"]] <- gMDP

            # show available classes in datasets
            TBStatus <- unique(colData(theObject)$TBStatus %>% as.character())

            # calculate the average gMDP score wihtin each class and sort from max to minimal
            gMDP_avg <- lapply(1:length(TBStatus), function(x){
              gMDP_class <- assays(theObject[,theObject$TBStatus == TBStatus[x]])[["gMDP"]]
              gMDP_class_avg <- apply(gMDP_class, 1, mean) %>% sort(decreasing = T)
              gMDP_class_avg
            })
            names(gMDP_avg) <- TBStatus

            # Use the original method in the paper
            # Identify genes with top 25% gMDP_avg within each class
            num_25 <- floor(nrow(gMDP)*0.25)
            gMDP_avg_gene <- lapply(gMDP_avg, function(x) x[1:num_25] %>% names())

            # claculate MDP of individual sample (sMDP)
            sMDP <- lapply(1:length(TBStatus), function(x,gMDP_avg_gene){
              gMDP_class <- assays(theObject[,theObject$TBStatus == TBStatus[x]])[["gMDP"]]
              sMDP_class <- gMDP_class[which(row.names(gMDP_class) %in% gMDP_avg_gene[[x]]),]
              sMDP_sample <- apply(sMDP_class, 2, mean)
              sMDP_sample
            }, gMDP_avg_gene)
            names(sMDP) <- TBStatus

            sMDP_data <- do.call(rbind, lapply(1:length(sMDP), function(x) {
              data.frame(Sample=names(sMDP[[x]]), sMDP = as.vector(sMDP[[x]]), TBStatus = names(sMDP[x]))
            }))

            return(sMDP_data)

          })

# Evaluate signatures using MDP score

setMethod("MDP", signature (theObject = "SummarizedExperiment", gset = "list"),
          function(theObject, gset = gset, assay_name = 1){
            # theObject data after normalization
            counts <- assays(theObject)[[assay_name]]
            # select the reference group, calculate mean and sd
            theObject_ref <- assays(theObject[,theObject$TBStatus == "Control"])[[assay_name]]
            mean_ref <- apply(theObject_ref, 1, mean)
            sd_ref <- apply(theObject_ref, 1, sd)

            # calculate gMDP for each gene i and sample s
            gMDP <- apply(counts, 2, function(x){(x-mean_ref)/sd_ref})

            # replace small gMDP values, set abs(gMDP) < 2 equals to 0
            gMDP[abs(gMDP)<2] <- 0

            # Include the gMDP assays in the existing SummarizedExperiment Object
            assays(theObject)[["gMDP"]] <- gMDP

            # show available classes in datasets
            # TBStatus <- unique(colData(theObject)$TBStatus %>% as.character())
            # Great methods, copy from gsva
            mapped.gset <- lapply(gset,function(x, y) na.omit(match(x, y)),
                                           row.names(gMDP))

            # claculate MDP of individual sample (sMDP)
            sMDP <- lapply(1:length(mapped.gset), function(x){
              # gMDP_class <- assays(theObject[,theObject$TBStatus == TBStatus[x]])[["gMDP"]]
              sMDP_class <- gMDP[mapped.gset[[x]],]
              sMDP_sample <- apply(sMDP_class, 2, mean)
              sMDP_sample
            })
            names(sMDP) <- names(mapped.gset)

            sMDP_sig <- do.call(cbind,sMDP) %>% data.frame()
            TBStatus <- colData(theObject)[,"TBStatus"]
            names(TBStatus) <- row.names(colData(theObject))

            sMDP_data <- cbind(sMDP_sig,TBStatus)
            return(sMDP_data)

          })

setMethod("MDP", signature (theObject = "MultiAssayExperiment", gset = "list"),
          function(theObject, gset = gset, experiment_type = "assay_reduce", assay_name = 1, GSE = GSE){

            # Create SummarizedExperiment
            col_data <-  colData(theObject)
            if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
              index <- sapply(1:length(colnames(theObject[[experiment_type]])), function (i)
                which(row.names(col_data) %in% colnames(theObject[[experiment_type]])[i]))

              col_data <- col_data[index,]
            }

            theObject <- SummarizedExperiment::SummarizedExperiment(assays = list(counts= as.matrix(theObject[[experiment_type]])),
                                                                    colData = col_data)

            # theObject data after normalization
            counts <- assays(theObject)[[assay_name]]
            # select the reference group, calculate mean and sd
            theObject_ref <- assays(theObject[,theObject$TBStatus == "Control"])[[assay_name]]
            mean_ref <- apply(theObject_ref, 1, mean) # 1 is running row by row
            sd_ref <- apply(theObject_ref, 1, sd)

            # calculate gMDP for each gene i and sample s
            # Get its absolute value
            gMDP <- apply(counts, 2, function(x){(x-mean_ref)/sd_ref}) # 2 is runnig column by column. Column vector calculation

            # replace small gMDP values, set abs(gMDP) < 2 equals to 0
            gMDP[abs(gMDP)<2] <- 0

            # set NaN value to 0, those with 0/0
            gMDP[is.nan(gMDP)] <- 0

            # set Inf value to 2, those with x/0
            gMDP[is.infinite(gMDP)] <- 2

            # Include the gMDP assays in the existing SummarizedExperiment Object
            assays(theObject)[["gMDP"]] <- gMDP

            # show available classes in datasets
            # TBStatus <- unique(colData(theObject)$TBStatus %>% as.character())
            # Great methods, copy from gsva
            mapped.gset <- lapply(gset,function(x, y) na.omit(match(x, y)),
                                  row.names(gMDP))

            ## remove gene sets from the analysis for which no features are available
            #gSetsLen <- sapply(mapped.gset,length)
            #gSet[gSetsLen >= 1 & gSetsLen <= Inf]

            # claculate MDP of individual sample (sMDP)
            sMDP <- lapply(1:length(mapped.gset), function(x){
              # gMDP_class <- assays(theObject[,theObject$TBStatus == TBStatus[x]])[["gMDP"]]
              sMDP_class <- gMDP[mapped.gset[[x]],]

              # Get mean value for each sample
              sMDP_sample <- apply(sMDP_class, 2, mean)
              sMDP_sample
            })
            names(sMDP) <- names(mapped.gset)

            sMDP_sig <- do.call(cbind,sMDP) %>% data.frame()
            TBStatus <- colData(theObject)[,"TBStatus"]
            names(TBStatus) <- row.names(colData(theObject))

            sMDP_data <- cbind(sMDP_sig,TBStatus, GSE)
            sMDP_data
            return(sMDP_data)

          })

