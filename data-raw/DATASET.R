## code to prepare `DATASET` dataset goes here

usethis::use_data("DATASET")

library(purrr)
library(devtools)
############ Function to add a list of summarized experiment data ############
add_new_data <- function(edit_files){
  library(purrr)
  library(devtools)
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

}

#----------------------------------------------------
# Add 3 data on Jan. 24th
edit_files <- paste0(c("GSE62147_sobject","GSE25534_sobject","GSE41055_sobject"),".RDS")
add_new_data(edit_files)

#----------------------------------------------------
# Modify RNA-seq data on Feb. 22th
add_new_data_mobject <- function(edit_files){
  library(purrr)
  library(devtools)
  sobject_total <- lapply(edit_files,
                          function(x) readRDS(paste0("~/Desktop/curatedTBData/data-raw/",x)))
  # sobject_total_name_ori <- lapply(sobject_total, function(x) names(assays(x)))
  # sobject_total_name <- lapply(sobject_total_name_ori, function(x) gsub("(.+?)(\\_.*)", "\\1", x)) %>% unlist()
  # sobject_total_name
  names(sobject_total) <- gsub("\\..*","",edit_files)

  purrr::walk2(sobject_total, paste0(names(sobject_total)), function(obj, name) {
    assign(name, obj)
    do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
  })

}

edit_files <- paste0(c("GSE101705_mobject", "GSE107104_mobject", "GSE89403_mobject",
                      "GSE79362_mobject","GSE94438_mobject", "GSE107991_mobject",
                      "GSE107992_mobject","GSE107993_mobject","GSE107994_mobject"),".RDS")
add_new_data_mobject(edit_files)

#----------------------------------------------------.
# Add 1 data on Mar. 1st
GSEXXXXX_sobject <- readRDS("data-raw/GSEXXXXX_sobject.RDS")
usethis::use_data(GSEXXXXX_sobject)

#----------------------------------------------------.
# Edit GSE112104 on Apr. 1st
add_new_data_mobject("GSE112104_mobject.RDS")

#----------------------------------------------------.
# Edit GSE107994,GSE79362 on Apr. 5th
add_new_data_mobject("GSE107994_mobject.RDS")
add_new_data_mobject("GSE79362_mobject.RDS")

#----------------------------------------------------.
# Edit GSE94438 on Apr. 5th
add_new_data_mobject("GSE94438_mobject.RDS")

#----------------------------------------------------.
# store Objects into indivual matrix

f1 <-  list.files("data-raw")
f2 <-  f1[grep("GSE",f1)]
total <- lapply(f2, function(x) readRDS(paste0("data-raw/",x)))

names(total) <- gsub("\\..*","",f2)

library(devtools)
purrr::walk2(total, names(total), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})

# Save Summary table
DataSummary <- readxl::read_excel("data-raw/Data_summaryforpackage.xlsx",
                                  sheet = "DataSummary")
SignatureInfo <- readxl::read_excel("data-raw/Data_summaryforpackage.xlsx",
                                  sheet = "SignatureInfo")
library(devtools)
use_data(DataSummary,compress = "xz", overwrite = TRUE)
use_data(SignatureInfo,compress = "xz", overwrite = TRUE)

load("~/Downloads/GSE79362_column_data.rda")
GSE79362_column_data$Age <- as.numeric(as.character(GSE79362_column_data$Age))
GSE79362_column_data$Gender <- ifelse(GSE79362_column_data$Gender=="female", "Female","Male")
GSE79362_column_data$Ethnicity <- as.character(GSE79362_column_data$Ethnicity)
GSE79362_column_data$Bin <- as.numeric(as.character(GSE79362_column_data$Bin))
GSE79362_column_data$TST <- as.numeric(as.character(GSE79362_column_data$TST ))
library(devtools)
use_data(GSE79362_column_data,compress = "xz", overwrite = TRUE)
saveRDS(GSE79362_column_data,"~/Desktop/curatedTBData/data-raw/GSE79362_column_data.RDS")


GSE79362_files <- list.files("data-raw",pattern = "GSE79362")
GSE79362_list <- lapply(GSE79362_files, function(x)
  readRDS(paste0("data-raw/",x)))

names(GSE79362_list) <- gsub("\\..*","",GSE79362_files)
library(devtools)
purrr::walk2(GSE79362_list, names(GSE79362_list), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})
