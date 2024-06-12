# install_packages.jl
import Pkg

# Activate the current project
Pkg.activate(".")  # or specify the path to your project

# List of packages to install
packages = [
    "JuMP", "HiGHS", "Ipopt",
    "LinearAlgebra", "DelimitedFiles", "Dates", "Logging",
    "BenchmarkTools",
    "Random", "Distributions",
    "BenchmarkProfiles",
    "Plots", "LaTeXStrings",
    "TimeZones",
    "DataFrames", "CSV", "XLSX",
    "SQLite"
]

# Install each package if not already installed
Pkg.add(url="https://github.com/mmogib/UIDFPA.jl.git")
for pkg in packages
    Pkg.add(pkg)
end

println("All packages are installed and ready to use.")