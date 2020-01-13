## code to prepare `DATASET` dataset goes here

usethis::use_data("DATASET")

library(purrr)
library(devtools)

sobject_files <- list.files(path = "data-raw",pattern = "*.RDS")

sobject_total <- lapply(sobject_files, function(x) readRDS(paste0("data-raw/",x)))
sobject_total_name_ori <- lapply(sobject_total, function(x) names(assays(x)))
sobject_total_name <- lapply(sobject_total_name_ori, function(x) gsub("(.+?)(\\_.*)", "\\1", x)) %>% unlist()
sobject_total_name
names(sobject_total) <- sobject_total_name

purrr::walk2(sobject_total, paste0(names(sobject_total), "_sobject"), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})

#############################################
edit_files <- paste0(c("GSE107991_sobject","GSE107992_sobject","GSE107993_sobject",
                       "GSE107994_sobject", "GSE62525_sobject", "GSE79362_sobject"),".RDS")
sobject_total <- lapply(edit_files,
                        function(x) readRDS(paste0("data-raw/",x)))
sobject_total_name_ori <- lapply(sobject_total, function(x) names(assays(x)))
sobject_total_name <- lapply(sobject_total_name_ori, function(x) gsub("(.+?)(\\_.*)", "\\1", x)) %>% unlist()
sobject_total_name
names(sobject_total) <- sobject_total_name

purrr::walk2(sobject_total, paste0(names(sobject_total), "_sobject"), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})
#Hiseq_sobject <- readRDS('data-raw/Hiseq_sobject.RDS')
#Illumina_sobject <- readRDS('data-raw/Illumina_sobject.RDS')
#Other_sobject <- readRDS('data-raw/Other_sobject.RDS')


#library(purrr)
#library(devtools)

#walk2(Other_sobject, names(Other_sobject), function(obj, name) {
#  assign(name, obj)
#  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
#})

#usethis::use_data(Affy_sobject)
#usethis::use_data(Hiseq_sobject)
#usethis::use_data(Illumina_sobject)
#usethis::use_data(Other_sobject)

