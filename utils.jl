function outputfolder(root::String="./results")
    zdt = now(tz"Asia/Riyadh")
    dayfolder = Dates.format(zdt, "yyyy_mm_dd")
    # hourfolder = Dates.format(zdt, "HH")
    root_dir = mkpath("$root/$dayfolder")
    return root_dir
end

function outputfilename(name::String, extension::String; root::String="./results", suffix::Union{Nothing,String}=nothing)
    root_dir = outputfolder(root)
    filename = if isnothing(suffix)
        "$root_dir/$name.$extension"
    else
        "$root_dir/$(name)_$suffix.$extension"
    end
    filename
end

function createStartingPoints()
    Dict(
        :onestimeshalf => n -> (1 / 2) * ones(n),
        :onestimes2 => n -> -2.0 * ones(n),
        :onestimes8 => n -> 8.0 * ones(n),
        :onestimes11halfs => n -> (11 / 2) * ones(n),
        :onestimes3halfs => n -> (3 / 2) * ones(n),
        :tenth => n -> 0.1 * ones(n),
        :negativetenth => n -> -0.1 * ones(n),
        :allones => n -> ones(n), # 8
        :negativeones => n -> -ones(n),
        :ntenthn => n -> begin
            i = Float64((n - 1)) / Float64(n)
            vcat(i, 0.1 * ones(n - 2), i)
        end,
        :fourn => n -> begin
            (1 / Float64((4.0 * Float64(n)^2))) * collect(1:n)
        end,
        :allzeros => n -> zeros(n),
        :tenzeros => n -> vcat(10, zeros(n - 1)),
        :ninezeros => n -> vcat(9, zeros(n - 1)),
        :threezeros => n -> vcat(3, zeros(n - 1)),
        :threeszeros => n -> [iseven(i) ? 0.0 : 3.0 for i in 1:n], #16
        :zerotwos => n -> vcat(0, 2 * ones(n - 1)),
        :zeroones => n -> vcat(0, ones(n - 1)),
        :oneszero => n -> vcat(ones(n - 1), 0),
        :halfn => n -> [1.0 / ((2.0)^i) for i in 1:n], #20
        :negativehalfn => n -> -[1.0 / ((2.0)^i) for i in 1:n],
        :nth => n -> [1 / i for i in 1:n],
        :negativenth => n -> -[1 / i for i in 1:n]
    )
end

function createPolyhedral(n::Int, m::Int=1; lower::Union{Number,Nothing}=nothing, upper::Union{Number,Nothing}=nothing, bvalue::Union{Number,Nothing}=nothing)
    A = ones(m, n)
    b = isnothing(bvalue) ? Float64.(n * ones(m)) : Float64.(bvalue * ones(m))
    lb = isnothing(lower) ? -Inf64 * ones(n) : lower * ones(n)
    ub = isnothing(upper) ? Inf64 * ones(n) : upper * ones(n)
    Polyhedral(A, b, lb, ub)
end




function plot_profiles(T, labels, xlabel)
    # T = hcat(T[:, 3], T[:, 6], T[:, 5], T[:, 4], T[:, 1], T[:, 2])
    plt_options = (line=(:dash, 2),
        size=(660, 450),
        palette=[:black, :blue, :green, :cyan, :red, :yellow], legend=:bottomright,
        frame=:box
    )
    # pythonplot()
    # ["SDM-AP", "SGM-AP1", "SGM-AP2", "SGM-AP3", "L-BFGS-AP", "MNM-AP"]
    p = performance_profile(PlotsBackend(), T, labels; xlabel, ylabel="Solved problems (%)", logscale=true, plt_options...)
    # fig_folder = outputfolder("paperCode/results")
    # savefig(p, "$fig_folder/gonacalev.png")
    return p
end


function solveWithNL(F::Function, x0::Vector{Number})
    sol = nlsolve(F, x0)
    sol.zero
end

function solveWithNL(problems, n::Int)
    F = problems[n][2].F
    points = map(last, problems[n][3])
    sols = map(x0 -> nlsolve(F, x0), points)
    sols
end

function createDB(; overwrite::Bool=false)

    db_name = "./results/results.db"
    result_tbl = "results"
    points_tbl = "initial_points_names"
    db = SQLite.DB(db_name)
    if overwrite
        SQLite.execute(db, "DROP TABLE IF EXISTS $result_tbl;")
    end
    SQLite.execute(
        db,
        """
CREATE TABLE IF NOT EXISTS $result_tbl (
    FunctionEvaluations INTEGER,
    FunctionNorm REAL,
    Iterations INTEGER,
    LSIterations INTEGER,
    Projections INTEGER,
    dim INTEGER,
    flag TEXT,
    line_search TEXT,
    message TEXT,
    problem TEXT,
    projection TEXT,
    search_direction TEXT,
    time REAL,
    x0 TEXT,
    insert_date TEXT,
    insert_time TEXT
);
"""
    )
    SQLite.execute(
        db,
        """
CREATE TABLE IF NOT EXISTS $points_tbl (
     id INTEGER PRIMARY KEY,
     name TEXT,
     value TEXT
);
"""
    )
    ipts_inserted = DBInterface.execute(db, "select * from $points_tbl;") |> DataFrame |> nrow
    if ipts_inserted == 0
        inps = getInitialPointsNames()
        ip_names = String.(keys(ips))
        ip_values = [values(ips)...]
        np_df = DataFrame(id=[i for i in 1:length(ip_names)], name=ip_names, value=ip_values)
        SQLite.load!(np_df, db, points_tbl)
    end
    DBInterface.execute(db, "CREATE VIEW IF NOT EXISTS results_v as SELECT t1.*,t2.value FROM $result_tbl t1 INNER JOIN $points_tbl t2 on t1.x0=t2.name;")
    db
end

function insertInDB(df::DataFrame; overwrite::Bool=false, result_table="results")
    db = createDB(overwrite=overwrite)
    SQLite.load!(df, db, result_table)
    SQLite.close(db)
end

function loadResultsFromDB()
    # insert_date::String, insert_time::String
    db = createDB()
    df = DBInterface.execute(db, "SELECT * FROM results_v;") |> DataFrame
    SQLite.close(db)
    df
end

function loadResultsFromDB(sqlwhere::String)
    # insert_date::String, insert_time::String
    db = createDB()
    df = DBInterface.execute(db, "SELECT * FROM results_v  where $sqlwhere") |> DataFrame
    SQLite.close(db)
    df
end

function deleteResultsFromDB(sqlwhere::String)
    # insert_date::String, insert_time::String
    df = loadResultsFromDB(sqlwhere::String)
    db = createDB()
    DBInterface.execute(db, "DELETE FROM results  where $sqlwhere")
    SQLite.close(db)
    df
end