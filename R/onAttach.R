".onAttach" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = pkg))
  title <- packageDescription(pkg, lib = mylib)$Title
  ver <- packageDescription(pkg, lib = mylib)$Version
  packageStartupMessage(paste(pkg, ": ", title, "\nVersion: ", ver, "\n", sep=""))
}
