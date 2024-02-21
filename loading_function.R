source("libraries.R")
list_file <- list.files(path = "functions/")

for (i in list_file){
  str_temp <- str_c("functions/",i)
  source(str_temp)
}
rm(list_file)
rm(str_temp)
rm(i)
