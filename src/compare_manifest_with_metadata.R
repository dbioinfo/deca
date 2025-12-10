library(tidyverse)

#first, navigate to the directory where the metadata file I sent you and the output of the ls command are
setwd("path/to/directory/")

#next read in each file, replace the text I have here with the names of the files
file1<-read_csv('PUT_METADATA_HERE.csv')
file2<-read_csv("PUT_LS_OUTPUT_HERE.csv", col_names = c("fname"))

#generate some stats
print(paste0("File Names in both lists: ", intersect(file1$fname, file2$fname)))
print(paste0("File Names exclusive to metadata: ", setdiff(file1$fname, file2$fname)))
print(paste0("File Names exclusive to file list: ", intersect(file1$fname, file2$fname)))
print(paste0("Proportion of files identified: ", intersect(file1$fname, file2$fname) / nrow(file2) ))