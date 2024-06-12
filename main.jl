include("init.jl")


"""
# 1) run experiments

"""
sols = run_experiement_1(
    problems_list=[i for i in 1:28],
    dim=[15_000, 50_000, 150_000],
    save_data=true,
    mxiters=2_000,
    time_limit=nothing,
    linesearchs=[LSI, LSII, LSIII, LSIV, LSV, LSVI, LSVII]
)


"""
# 2) load data from database

"""
results_df = loadResultsFromDB()

"""
# 3) Make Plots

"""
plts = plotData()

"""
# 5) Failed Problems

"""
plts = plotFailedProblems(; save_plot=true, dim=100_00)
