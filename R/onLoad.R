".onLoad" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = "eco"))
  title <- packageDescription("eco", lib = mylib)$Title
  ver <- packageDescription("eco", lib = mylib)$Version
  cat(title, "\nVersion", ver, "\n")
}

