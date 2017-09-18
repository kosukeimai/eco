".onAttach" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = pkg))
  title <- packageDescription(pkg)$Title
  ver <- packageDescription(pkg)$Version
  packageStartupMessage(paste(pkg, ": ", title, "\nVersion: ", ver, "\n", sep=""))
}
