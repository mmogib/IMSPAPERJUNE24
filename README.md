# Paper Experiments
> Unified Inertial Derivative-Free Projection Methods for Constrained Nonlinear Equations 
>
> Authors:
> - Abdulkarim Hassan Ibrahim, 
> - Mohammed Alshahrani, 
> - Suliman Al-Homidan 

This repository contains the code to reproduce the results for the `UIDFPAF` algorithm on 28 test problems. Follow the steps below to set up the environment and run the experiments.

## Step 1: Install Julia

1. Download and install Julia from the [official website](https://julialang.org/downloads/).
2. Add Julia to your system PATH if it is not done automatically during installation.

## Step 2: Install Required Packages

1. Clone this repository or download the files to your local machine.
2. Open a terminal and navigate to the directory containing the repository.
3. Create a new Julia environment or activate an existing one by running:
    ```sh
    julia --project=.
    ```
5. Exit from Julia `REPL` by running
6. ```julia
   exit()
   ```
7. Run the package installation script to install all required packages:
    ```sh
    julia install_packages.jl
    ```

## Step 3: Run the Code

1. Run the main script to perform the experiments and generate results:
    ```sh
    julia main.jl
    ```

## Content of `main.jl`

The `main.jl` file is the main script that orchestrates the experiments, data loading, and plotting. Below is a detailed explanation of its contents:

### 1. Initialize

The script begins by including the initialization file:
```julia
include("init.jl")
```

### 2. Run Experiments

The `run_experiment_1` function runs the experiments on the specified problems and dimensions, saves the data, and sets the maximum iterations and time limits (nothing means no time limit).
```julia
sols = run_experiement_1(
    problems_list=[i for i in 1:28],
    dim=[15_000, 50_000,150_000],
    save_data=true,
    mxiters=2_000,
    time_limit=nothing,
    linesearchs=[LSI, LSII, LSIII, LSIV, LSV, LSVI, LSVII]
)
```

### 3. Load Data from Database

The `loadResultsFromDB` function loads the results data from the database (saved in `results` folder) into a DataFrame for analysis.
```julia
results_df = loadResultsFromDB()
```

### 4. Make Plots

The `plotData` function generates plots based on the results data.
```julia
plts = plotData()
```

### 5. Plot Failed Problems

The `plotFailedProblems` function generates plots for the problems that failed during the experiments.
```julia
plts = plotFailedProblems(; save_plot=true, dim=1000)
```

## Additional Information

For more details on the functions and their implementations, refer to the respective `.jl` files in the repository.
```

This `README.md` file provides a clear and detailed guide for anyone who wants to reproduce the results, including installing Julia, setting up the environment, running the experiments, and understanding the main script's content.