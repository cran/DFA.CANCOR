
.onAttach<- function(libname, pkgname){
  packageStartupMessage("**************************************************************************************************\n",
                        pkgname," ",packageDescription("DFA.CANCOR")$Version,
                        "\n\nPlease contact Brian O'Connor at brian.oconnor@ubc.ca if you have questions or suggestions.\n",
                        "**************************************************************************************************", 
                        appendLF = TRUE)
}

NULL


# cat("Package: DFA.CANCOR    Version: ", packageDescription("DFA.CANCOR")$Version,"\n")


# # .onAttach<- function(libname, pkgname){
  # packageStartupMessage("**************************************************************************************************\nWelcome to ",
                        # pkgname," ",packageDescription("DFA.CANCOR")$Version,
                        # "\n\nPlease contact Brian O'Connor at brian.oconnor@ubc.ca if you have questions or suggestions.\n",
                        # "**************************************************************************************************", 
                        # appendLF = TRUE)
# }

