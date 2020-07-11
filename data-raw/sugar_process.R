## code to prepare `sugar-process` dataset goes here

sugar_process <- read.csv("ash_color.txt", header=TRUE)

excitation_230 <- read.csv("excitation_230.txt", header=TRUE, check.names = FALSE)
excitation_240 <- read.csv("excitation_240.txt", header=TRUE, check.names = FALSE)
excitation_255 <- read.csv("excitation_255.txt", header=TRUE, check.names = FALSE)
excitation_290 <- read.csv("excitation_290.txt", header=TRUE, check.names = FALSE)
excitation_305 <- read.csv("excitation_305.txt", header=TRUE, check.names = FALSE)
excitation_325 <- read.csv("excitation_325.txt", header=TRUE, check.names = FALSE)
excitation_340 <- read.csv("excitation_340.txt", header=TRUE, check.names = FALSE)

sugar_process$excitation_230 <- as.matrix(excitation_230)
sugar_process$excitation_240 <- as.matrix(excitation_240)
sugar_process$excitation_255 <- as.matrix(excitation_255)
sugar_process$excitation_290 <- as.matrix(excitation_290)
sugar_process$excitation_305 <- as.matrix(excitation_305)
sugar_process$excitation_325 <- as.matrix(excitation_325)
sugar_process$excitation_340 <- as.matrix(excitation_340)


usethis::use_data(sugar_process, overwrite = TRUE)
