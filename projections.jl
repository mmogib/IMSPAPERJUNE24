
function projectOnBoxCheck(x::Vector{<:Real}; bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    if isnothing(bounds)
        true
    else
        l, u = bounds
        all(l .<= x .<= u)
    end
end
function projectOnBox(x::Vector{<:Real}; bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    if isnothing(bounds)
        x
    else
        l, u = bounds
        x .|> t -> min(t, u) .|> t -> max(t, l)
    end
end

function projectOnHalfSpaceCheck(x; β::Union{Nothing,Real}=nothing)
    n = length(x)
    b = isnothing(β) ? n : β
    sum(x) <= b
end

function projectOnHalfSpace(x; β::Union{Nothing,Real}=nothing)
    # % project an n-dimensional vector x to the the halfspace Sᵦ
    # % Sᵦ = { x : x ∈ Rⁿ, sum(x) <= β}
    # % Coded by
    n = length(x)
    b = isnothing(β) ? n : β
    if (sum(x) <= b)
        return x
    end

    bget = false
    s = sort(x, rev=true)
    tmpsum = 0
    for ii = 1:n-1
        tmpsum = tmpsum + s[ii]
        tmax = (tmpsum - b) / ii
        if tmax >= s[ii+1]
            bget = true
            break
        end
    end
    tmpsum = tmpsum + s[n]
    if ~bget
        tmax = (tmpsum - b) / n
    end

    P = if (tmpsum <= b)
        x
    else
        x .- tmax .|> r -> max(r, 0)
    end
    P
end
function projectOnTriangleCheck(x; β::Union{Nothing,Real}=nothing, lower=0)
    projectOnBoxCheck(x, bounds=(lower, Inf64)) && projectOnHalfSpaceCheck(x, β=β)
end
function projectOnTriangle(x; β::Union{Nothing,Real}=nothing, lower=0)
    # % project an n-dim vector x to the convex set Cn
    # % Cn = { x : x n-dim, x >= lower, sum(x) <= β}
    P = projectOnBox(x, bounds=(lower, Inf64)) |> p -> projectOnHalfSpace(p, β=β)
    P
end
function projectOnTriangleCheck(x; β::Union{Nothing,Real}=nothing, bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    projectOnBoxCheck(x, bounds=bounds) && projectOnHalfSpaceCheck(x, β=β)
end
function projectOnTriangle(x; β::Union{Nothing,Real}=nothing, bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    # % project an n-dim vector x to the convex set Cn
    # % Cn = { x : x n-dim, x >= lower, sum(x) <= β}
    P = projectOnBox(x, bounds=bounds) |> p -> projectOnHalfSpace(p, β=β)
    P
end

function projectOnRn(x)
    x
end