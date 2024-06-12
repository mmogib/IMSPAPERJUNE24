function getLineSearchNames()
    Dict(
        "LSI" => LSI,
        "LSII" => LSII,
        "LSIII" => LSIII,
        "LSIV" => LSIV,
        "LSV" => LSV,
        "LSVII" => LSVII,
        "LSVI" => LSVI,
    )
end
function getInitialPointsNames()
    initial_points_names = Dict(
        :tenth => "x₁",
        :negativetenth => "x₂",
        :allones => "x₃",
        :negativeones => "x₄",
        :ntenthn => "x₅",
        :allzeros => "x₆",
        :tenzeros => "x₇",
        :ninezeros => "x₈",
        :threezeros => "x₉",
        :threeszeros => "x₁₀",
        :zerotwos => "x₁₁",
        :zeroones => "x₁₂",
        :oneszero => "x₁₃",
        :halfn => "x₁₄",
        :negativehalfn => "x₁₅",
        :nth => "x₁₆",
        :negativenth => "x₁₇"
    )
    initial_points_names
end
function plotFailedProblems(; save_plot::Bool=true, dim::Int=100_000)
    file_name = "failed_problems"
    plts = plotFailedProblems(dim)
    if save_plot
        for plt in plts
            for (name, fig) in plt
                output_file_svg = outputfilename("$(file_name)_$name", "svg")
                output_file_png = outputfilename("$(file_name)_$name", "png")
                Plots.svg(fig, output_file_svg)
                Plots.png(fig, output_file_png)
            end
        end
    end
    plts
end
function plotFailedProblems(dim::Int=100_000)
    initial_points_names = getInitialPointsNames()
    pnorms = gathernorms(dim=dim)
    mcolors = [:dimgray, :bisque, :coral, :darkolivegreen]
    map(pnorms) do p
        pnam = p[:problem_name]
        map(p[:ls_norms]) do ls
            lsname = ls[:ls_name]
            xsnorms = filter(x -> x[:solved] == false, ls[:xs_norms])
            xs = map(enumerate(xsnorms)) do (i, xall)
                xi = xall[:norms]
                xiname = String(xall[:x0_name])
                lstvalue = round(xi[end], digits=6)
                (xiname, scatter(xi, markersize=2, markercolor=mcolors[i], annotations=[(length(xi) / 2, 0.9 * maximum(xi), L"\psi(x_{%$(length(xi))})=%$lstvalue", 10)]))
            end
            # println(length(xs), " :: ", pnam * " - " * lsname)
            # map(x -> x[2], xs)
            if length(xs) == 0
                nothing
            else
                xslen = length(xs)
                xsw = Int(floor(xslen / 2))
                xsl = Int(ceil(xslen / 2))
                lbls = reshape(map(x -> initial_points_names[Symbol(x[1])], xs), 1, xslen)
                ylbls = L"\|\psi(x_k)(\|"
                xlbl = L"k"
                (pnam * " - " * lsname,
                    plot(map(x -> x[2], xs)..., layout=(2, 2), xlabel=xlbl, ylabel=ylbls, title=lbls, legend=:none))
            end
            # T = stack(xs)
            # Tname = map(x -> x[:x0_name], xsnorms)
        end |> filter(x -> !isnothing(x))
    end
end
function gathernorms(; dim::Int=100_000)
    # initial_points_names = getInitialPointsNames()
    # line_search_names = getLineSearchNames()
    sql = "SELECT DISTINCT problem from results  WHERE flag not in ('error','solved') order by length(problem),problem"
    db = createDB()
    data = DBInterface.execute(db, sql) |> DataFrame
    SQLite.close(db)
    problems = init_experiment(dim) |> filter(x -> x[1] in data[!, 1]) |> y -> map(x -> (x[1], x[2], x[3]), y)
    gathernorms.(problems)
end

function gathernorms(problem::Tuple{String,NECProblem,Vector{Tuple{Symbol,Vector{Float64}}}})
    (name, p, xs) = problem
    linesearchs = [LSI, LSII, LSIII, LSIV, LSV, LSVI, LSVII]
    params = UIDFPAParams(0.6, 0.25, 0.01, (x = 0.5) -> x)
    normsForProblem = map(linesearchs) do ls
        lsname = String(Symbol(ls))
        options = UIDFPAOptions(params, ls, compute_search_direction, 1e-6, 2000, 180, true, false, true)
        printstyled("Solving problem $name with lines search $lsname and search_direction compute_search_direction and $(length(xs)) points of size $dim ...\n", color=:blue, bold=true)
        @debug "Solving.. problem $name with lines search $lsname and $(length(xs)) points of size $n ...\n"
        normsforXs = map(xs) do (x0name, x0)
            x0dim = length(x0)
            printstyled("running $x0name of size $(x0dim).\n", color=:blue)
            sol = uidfpa(p, x0; options)
            norms = filter(x -> x != -1, sol.solution.FunctionNorms)
            Dict(
                :norms => norms,
                :x0_name => x0name,
                :dim => x0dim,
                :x0 => x0,
                :fx0norm => norm(p.F(x0)),
                :solved => isa(sol, SuccessSolution) ? true : false
            )
        end
        Dict(
            :ls_name => lsname,
            :xs_norms => normsforXs
        )
    end
    Dict(
        :problem_name => name,
        :ls_norms => normsForProblem
    )

end


function run_algorithm(n::Union{Int,Vector{Int}};
    linesearchs::Union{Function,Vector{Function}}=[LSIII, LSIV, LSVI],
    search_directions::Union{Function,Vector{Function}}=compute_search_direction,
    ϵ::Float64=1e-6, mxitrs::Int=1000,
    params::UIDFPAParams=UIDFPAParams(0.5, 0.25), choose_problems::Function=(i -> i > 5 || i < 3),
    is_inertia::Bool=true,
    is_approximate::Bool=true,
    time_limit::Union{Nothing,Int}=nothing
)
    headers = [:line_search, :search_direction, :x0, :problem, :flag, :message, :FunctionNorm, :Iterations, :Projections, :LSIterations, :FunctionEvaluations, :time, :dim, :projection]
    problem_sizes = isa(n, Int) ? [n] : n
    lss = isa(linesearchs, Function) ? [linesearchs] : linesearchs
    sdir = isa(search_directions, Function) ? [search_directions] : search_directions
    approximate = is_approximate ? :approximate : :exact
    gather_norms, tm_limit = isnothing(time_limit) ? (false, 0) : (true, time_limit)
    dfsall = map(problem_sizes) do n
        problems = init_experiment(n)

        dfs = map(lss) do ls
            lsname = String(Symbol(ls))
            dfsd = map(sdir) do sd
                sdname = String(Symbol(sd))
                options = UIDFPAOptions(params, ls, sd, ϵ, mxitrs, tm_limit, is_inertia, is_approximate, gather_norms)
                df1 = map(problems[filter(choose_problems, 1:length(problems))]) do (name, p, xs)
                    printstyled("Solving problem $name with lines search $lsname and search_direction $sdname and $(length(xs)) points of size $n ...\n", color=:blue, bold=true)
                    @debug "Solving.. problem $name with lines search $lsname and $(length(xs)) points of size $n ...\n"
                    df2 = map(xs) do (x0name, x0)
                        x0dim = length(x0)
                        printstyled("running $x0name of size $(x0dim).\n", color=:blue)
                        T = Matrix{Union{String,Symbol,<:Number,Nothing}}(undef, 10, 1)
                        col = 1
                        try

                            timed_sol = @timed uidfpa(p, x0; options)

                            sol = timed_sol.value
                            if isa(sol, SuccessSolution)
                                # println("solution", sol.solution.x)
                                T[col] = "solved" # flag
                                T[col+1] = sol.message
                                T[col+2] = sol.solution.Fnorm
                                T[col+3] = sol.solution.Iterations
                                T[col+4] = sol.solution.ProjectionIterations
                                T[col+5] = sol.solution.LineseatchIterations
                                T[col+6] = sol.solution.FunctionEvaluations
                                T[col+7] = timed_sol.time
                                T[col+8] = x0dim
                            else
                                T[col] = "failed" # flag
                                T[col+1] = sol.message
                                T[col+2] = sol.solution.Fnorm
                                T[col+3] = sol.solution.Iterations
                                T[col+4] = sol.solution.ProjectionIterations
                                T[col+5] = sol.solution.LineseatchIterations
                                T[col+6] = sol.solution.FunctionEvaluations
                                T[col+7] = timed_sol.time
                                T[col+8] = x0dim
                            end
                        catch er
                            T[col] = "error" # flag
                            T[col+1] = String(Symbol(er))
                            T[col+2] = NaN
                            T[col+3] = NaN
                            T[col+4] = NaN
                            T[col+5] = NaN
                            T[col+6] = NaN
                            T[col+7] = NaN
                            T[col+8] = x0dim
                        finally
                            T[col+9] = String(approximate)
                        end
                        DataFrame(Dict(zip(headers, vcat(lsname, sdname, String(Symbol(x0name)), name, T))))
                    end
                    vcat(df2...)
                end
                vcat(df1...)
            end
            vcat(dfsd...)
        end
        vcat(dfs...)
    end
    return vcat(dfsall...)
end

# 1. save data
function save_experiments(df::DataFrame; excel=false, file_name::String="experiment_1", suffix::Union{String,Nothing}=nothing)
    if excel
        output_file_csv = outputfilename(file_name, "xlsx"; suffix)
        XLSX.writetable(output_file_csv, df, overwrite=true, sheetname=file_name)
    else
        output_file_csv = outputfilename(file_name, "csv"; suffix)
        CSV.write(output_file_csv, df)
    end
end





# 2. create profiles
function save_plots(df::DataFrame; linesearchs::Vector{Symbol}=[:LSIII, :LSIV, :LSVI], file_name::String="experiment_1", suffix::Union{String,Nothing}=nothing)
    colors = [:red, :darkorange1, :green, :blue, :purple, :black, :yellow]
    (w, h) = Plots._plot_defaults[:size]
    ys = 0:0.1:1.01

    plts = map([(:Iterations, "Iterations"), (:time, "Time"), (:FunctionEvaluations, "Function Evaluations")]) do (data_point, data_header)
        T = Matrix{Float64}(undef, Int(DataFrames.nrow(df) / length(linesearchs)), length(linesearchs))
        foreach(enumerate(linesearchs)) do (i, name)
            d = select(filter(r -> r[:line_search] == String(name), df), [data_point, :flag] => ByRow((itr, f) -> f == "solved" ? itr : NaN) => Symbol(data_header))
            T[:, i] = d[!, Symbol(data_header)]
        end
        pdata = performance_profile_data(T)
        max_ratio = pdata[3]
        xs = 1:ceil(max_ratio + 0.5)
        p = performance_profile(PlotsBackend(), T, String.(linesearchs);
            title="$(data_header)",
            logscale=true,
            size=(1.2w, h),
            xlabel="Performance ratio",
            ylabel="Solved problems (%)",
            legendfontsize=8,
            linestyle=:dash,
            palette=colors,
            linewidth=2.5,
            minorgrid=true, leg=:bottomright
        )
        plot(p, xticks=(xs, map(x -> "$(Int(x))", xs)),
            yticks=(ys, map(x -> "$x", ys))), data_header
    end

    for (p, name) in plts
        output_file_svg = outputfilename("$(file_name)_$name", "svg"; suffix)
        output_file_png = outputfilename("$(file_name)_$name", "png"; suffix)

        Plots.svg(p, output_file_svg)
        Plots.png(p, output_file_png)
    end
end

function run_experiement_1(n::Union{Int,Vector{Int}};
    linesearchs::Union{Function,Vector{Function}}=[LSIII, LSIV],
    search_directions::Union{Function,Vector{Function}}=compute_search_direction,
    ϵ::Float64=1e-6,
    mxitrs::Int=1000,
    params::UIDFPAParams=UIDFPAParams(0.5, 0.25),
    saveit::Bool=true,
    plotit::Bool=true,
    choose_problems::Function=(i -> i > 5 || i < 3),
    is_inertia::Bool=true,
    is_approximate::Bool=true,
    problemsRange::Union{Vector{Int},UnitRange}=1:18,
    time_limit::Union{Nothing,Int}=nothing
)
    # 1 solve
    printstyled("Starting \n", color=:blue, bold=true)
    df = run_algorithm(n; linesearchs, search_directions, ϵ, mxitrs, params, choose_problems, is_inertia, is_approximate, time_limit)
    problems_str = join(filter(i -> choose_problems(i), problemsRange) .|> j -> "p_$j", "_")
    suffix = isa(linesearchs, Vector) ? join(linesearchs, "_") * problems_str : String(Symbol(linesearchs)) * "_" * problems_str
    # 2 save
    if (saveit)
        printstyled("==========================================================================\n", color=:blue, bold=true)
        printstyled("Saving data ...... \n", color=:blue, bold=true)
        # linesearchs = isa(linesearchs, Function) ? [linesearchs] : linesearchs
        save_experiments(df; suffix)
    end
    # 3 plot and save
    if plotit
        printstyled("==========================================================================\n", color=:blue, bold=true)
        printstyled("Plotting data ...... \n", color=:blue, bold=true)
        save_plots(df; linesearchs=Symbol.(linesearchs), suffix)
    end
    zdt = now(tz"Asia/Riyadh")
    insdate = Dates.format(zdt, "yyyy_mm_dd")
    instime = Dates.format(zdt, "HH_MM")
    hcat(df, DataFrame(insert_date=repeat([insdate], nrow(df)), insert_time=repeat([instime], nrow(df))))
end

function runExperiement1(linesearchs, problemsRange, search_directions=compute_search_direction; dim=[1000], mxiters::Int=2000,
    time_limit::Union{Nothing,Int}=180)
    debuglogger = ConsoleLogger(stderr, Logging.Info)
    df = DataFrame()
    sols = with_logger(debuglogger) do
        for j in problemsRange
            printstyled("Solving problem $j \n", color=:blue, bold=true)
            printstyled("==========================================================================\n", color=:voilet, bold=true)
            df0 = run_experiement_1(dim;
                linesearchs,
                search_directions,
                ϵ=1e-6,
                mxitrs=mxiters,
                params=UIDFPAParams(0.6, 0.25, 0.01, (x = 0.5) -> x),
                saveit=false,
                plotit=false,
                choose_problems=(i -> i == j),
                is_inertia=true,
                is_approximate=false,
                problemsRange=problemsRange,
                time_limit=time_limit
            )
            df = vcat(df, df0)
        end
        df
    end
    sols
end
nan2empty(x) = isnan(x) ? "" : x
nan2empty(x::Missing) = ""
function run_experiement_1(; problems_list::Vector{Int}=collect(1:28),
    dim::Vector{Int}=[15_000, 50_000, 150_000],
    save_data::Bool=true,
    save_for_ibrahim::Bool=true,
    mxiters::Int=2000,
    time_limit::Union{Nothing,Int}=180,
    linesearchs=[LSI, LSII, LSIII, LSIV, LSV, LSVI, LSVII]
)
    ## DB SETTING
    db_name = "./results/paper3.db"
    db_table = "results"
    # insert_date, insert_time = "2024_05_30", "14_01"
    # overwrite = false


    # linesearchs = LSV
    # search_directions = [steepest_descent, cruz_raydan, ye_zhou, abubaker_mohammad, compute_search_direction]


    # disable_logging(Logging.Info)

    sols = runExperiement1(linesearchs, problems_list; dim=dim, mxiters=mxiters, time_limit=time_limit)
    deleteddf = DataFrame()
    plist = map(x -> "'p$x'", problems_list) |> x -> join(x, ",")
    sqlwhere = "problem in ($plist)"
    if save_data
        deleteddf = deleteResultsFromDB(sqlwhere)
        insertInDB(sols; overwrite=false)
    end
    ibrahim_df = DataFrame()
    if save_for_ibrahim
        linesearchs_str = String.(Symbol.(linesearchs))

        stored_dfs = loadResultsFromDB(sqlwhere)

        df_colums = map(enumerate(linesearchs_str)) do (i, lsstr)
            # (lsstr, "dim", "x0", "iter", "Feval", "time")
            df = stored_dfs[stored_dfs.line_search.==lsstr, :]
            DataFrame(
                lsstr * "_dim" => df[!, :dim],
                lsstr * "_problem" => df[!, :problem],
                lsstr * "_x0" => df[!, :x0],
                lsstr * "_Fnorm" => map(nan2empty, df[!, :FunctionNorm]),
                lsstr * "_iter" => map(nan2empty, df[!, :Iterations]),
                lsstr * "_Feval" => map(nan2empty, df[!, :FunctionEvaluations]),
                lsstr * "_time" => map(nan2empty, df[!, :time])
            )
        end |> df -> hcat(df...)
        save_experiments(df_colums; excel=true, file_name="experiment_1_problem_1_32", suffix="exact_projection_with_inertia")
        ibrahim_df = df_colums
    end


    Dict(:new_df => sols, :old_df => deleteddf, :ibrahim_df => ibrahim_df)

end

function plotData()
    println("Plotting...")
    linesearchs = [LSI, LSII, LSIII, LSIV, LSV, LSVI, LSVII]
    stored_dfs = loadResultsFromDB()
    save_plots(stored_dfs; linesearchs=Symbol.(linesearchs))
end