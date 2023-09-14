source("libraries.R")

list_file <- list.files(path = "fun_confirmed/")

for (i in list_file){
  str_temp <- str_c("fun_confirmed/",i)
  source(str_temp)
}
rm(list_file)
rm(str_temp)
rm(i)
