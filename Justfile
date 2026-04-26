# Local Variables:
# mode: makefile
# End:

set shell := ["bash", "-cu"]

alias fc := fmt-check
alias t := test
alias f := fmt

# list commands
default:
	just --list

# format with Runic
fmt:
	julia  --startup-file=no --project=@runic -m Runic --inplace src/

# check format with Runic
fmt-check:
	julia --startup-file=no --project=@runic -m Runic --check --diff src/

# run Pkg.test
test:
	julia --project=. --startup-file no -e "import Pkg; Pkg.test()"

# install Runic in @runic
install-runic:
	julia  --startup-file=no --project=@runic -e 'using Pkg; Pkg.add("Runic")'
