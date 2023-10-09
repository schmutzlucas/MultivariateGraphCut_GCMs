#' Install and Load Necessary Libraries
#'
#' This function reads a list of required R packages from a text file and installs any packages that are not already installed. It then loads these packages into the R environment.
#'
#' @param package_list_file A character string specifying the path to the text file that contains the list of required packages. The default is "package_list.txt".
#'
#' @return NULL. The function loads the necessary libraries into the R environment.
#'
#' @examples
#' # Example usage:
#' # install_and_load_libraries("package_list.txt")
#'
#' @import utils
#' @export
install_and_load_libraries <- function(package_list_file = "package_list.txt") {
  list_of_packages <- read.table(package_list_file, sep = "\n")$V1
  new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]

  if (length(new.packages)) {
    install.packages(new.packages, repos = "https://cloud.r-project.org")
  }

  lapply(list_of_packages, library, character.only = TRUE)
}
