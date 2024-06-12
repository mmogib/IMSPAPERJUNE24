"""
Source:
Problem 4.6 in
Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, 
Modified optimal Perry conjugate gradient method for solving system of monotone equations with applications,Applied Numerical Mathematics,
Volume 184, 2023, Pages 431-445, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2022.10.016.

Projected on [0, ∞] or [-1, ∞]

"""
function OPolynomialSineCosine(x)
    N = 1 / length(x)
    s = N * sum(x)
    x .* cos.(x .- N) .* (
        sin.(x) .- 1 .- (x .- 1) .^ 2 .- s
    )
end

"""
Source: 
Problem 5.2
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96. https://www.global-sci.org/intro/article_detail/jcm/8675.html#
Example: 4.1 in
Wang, C., Wang, Y. & Xu, C. A projection method for a system of nonlinear monotone equations with convex constraints. 
Math Meth Oper Res 66, 33–46 (2007). https://doi.org/10.1007/s00186-006-0140-y
"""
function OExponetialI(x)
    exp(x) - 1
end



"""
Source:
Problem 3 in Appendix
William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

"""
function OExponetialIII(x)
    X2 = x .^ 2
    Ii = [i / 10 for i in 1:length(x)]
    X21NED = X2[1:end-1] + exp.(-X2[1:end-1])
    Y = vcat(X21NED, exp(-X2[end]))
    Ii .* (1 .- Y)
end


"""
Source:
Problem 9 
J. Sabi’u, A. Shah, M.Y. Waziri, M.K. Dauda, A new hybrid approach for solving large-scale monotone nonlinear equations, 
J. Math. Fund. Sci. 52 (2020) 17–26 https://doi.org/10.5614/j.math.fund.sci.2020.52.1.2


"""
function OPolynomialI(x)
    vcat(
        4 .* x[1:end-1] .+ (x[2:end] .- 2 .* x[1:end-1]) .- ((x[2:end] .^ 2) / 3),
        4 * x[end] + (x[end-1] - 2 * x[end]) - x[end-1]^2 / 3
    )
end



"""
Source:
Problem 2 in
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.

Project to Rn₊

Modified Problem 1. in
W.J. Zhou, D.H. Li, A globally convergent BFGS method for nonlinear monotone equations without any merit functionsO, 
Math. Comput. 77 (264) (2008) 2231–2240.
Projected on [-2, ∞]


"""
function OSmoothSine(x)
    2 * x - sin(x)
end



"""
Source:
Problem 1 in
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.tional and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

Project to Rn₊
"""
function ONonsmoothSine(x)
    2 * x - sin(abs(x))
end


"""
Source:
Problem 2 in
Zhensheng Yu, Ji Lin, Jing Sun, Yunhai Xiao, Liying Liu, Zhanhui Li, Spectral gradient projection method for monotone nonlinear equations with convex constraints,
Applied Numerical Mathematics, Volume 59, Issue 10, 2009, Pages 2416-2423, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2009.04.004.

with Projection on 
    C = {x ∈ Rⁿ | ∑xᵢ ≤ n, xᵢ≥ -1}

Modified from Problem 1 in
Li Zhang, Weijun Zhou, Spectral gradient projection method for solving nonlinear monotone equations, 
Journal of Computational and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

"""
function OModifiedNonsmoothSine(x)
    x - sin(abs(x - 1))
end

"""
Source:
Problem 2
Gao PT, He CJ (2018) An efficient three-term conjugate gradient method for nonlinear monotone equations
with convex constraints. Calcolo. https://doi.org/10.1007/s10092-018-0291-2


with Projection on 
    C = {x ∈ Rⁿ | ∑xᵢ ≤ n, xᵢ≥ -1}

Modified from Problem 1 in
Li Zhang, Weijun Zhou, Spectral gradient projection method for solving nonlinear monotone equations, 
Journal of Computational and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

"""
function OModifiedNonsmoothSine2(x)
    x - sin(abs(x) - 1)
end


"""
Source:
Problem 1
Gao, P., He, C. An efficient three-term conjugate gradient method for nonlinear monotone equations with convex constraints. 
Calcolo 55, 53 (2018). https://doi.org/10.1007/s10092-018-0291-2

with Projection on 
    Rn+

"""
function OExponetialSineCosine(x)
    exp(2x) + 3 * sin(x) * cos(x) - 1
end

"""
Source: 
Problem 4.2 in
Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427,
https://doi.org/10.1016/j.cam.2020.112781.
Modified from Problem 3
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.
"""
function OModifiedTrigI(x)
    vcat(
        x[1] + sin(x[1]) - 1,
        -x[1:end-2] + 2 * x[2:end-1] + sin.(x[2:end-1]) .- 1,
        x[end] + sin(x[end]) - 1
    )
end

"""
Source:
Problem 10 in 
Y. Bing, G. Lin, An efficient implementation of Merrill’s method for sparse or partially separable systems of nonlinear equations, 
SIAM. J. Optim.1 (2) (1991) 206–221.
"""
function OTridiagonal(x)
    Ii = [i for i in 1:length(x)]
    x0 = vcat(0, x)
    x1 = vcat(x, 0)
    x .- exp.(cos.(Ii .* (x0[1:end-1] .+ x .+ x1[2:end])))
end


"""
Source:
Problem 4.4 in
Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427,
https://doi.org/10.1016/j.cam.2020.112781.
Modified from Problem 10 in 
Y. Bing, G. Lin, An efficient implementation of Merrill’s method for sparse or partially separable systems of nonlinear equations, 
SIAM. J. Optim.1 (2) (1991) 206–221.
"""
function OModifiedTridiagonal(x)
    n = length(x)
    vcat(
        x[1] - exp(cos((x[1] + x[2]) / (n + 1))),
        -x[2:end-1] - exp.(cos.((x[1:end-2] .+ x[2:end-1] .+ x[3:end]) ./ (n + 1))),
        x[end] - exp(cos((x[end-1] + x[end]) / (n + 1)))
    )
end
function OModifiedTridiagonalV2(x)
    n = length(x)
    F = zeros(n)

    F[1] = x[1] - exp(cos((x[1] + x[2]) / (n + 1)))

    for i in 2:n-1
        F[i] = -x[i] - exp(cos((x[i-1] + x[i] + x[i+1]) / (n + 1)))
    end

    F[n] = x[n] - exp(cos((x[n-1] + x[n]) / (n + 1)))

    return F
end

"""
Source:
Problem 10 in Appendix
William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493


Problem 4.5 in
Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2020.112781

Projected on [-1, ∞]

"""
function OLogarithmic(x)
    log.(x .+ 1) .- (x / length(x))
end


"""
Source:
Problem 10 in Appendix
William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493


Modified in Problem 4.2 in
Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, 
Modified optimal Perry conjugate gradient method for solving system of monotone equations with applications,Applied Numerical Mathematics,
Volume 184, 2023, Pages 431-445, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2022.10.016.

Projected on [0, ∞] or [-1, ∞]

"""
function ONonmoothLogarithmic(x)
    log.(abs.(x) .+ 1) .- (x / length(x))
end

function Oding2017FunII(x)
    x - sin(abs(x - 1))
end

function Oding2017Diagonal6Fun(x)
    cp5(x)
end

# see https://www.cuter.rl.ac.uk//Problems/classification.shtml for classification

# NAME: ARWHEAD         see (https://bitbucket.org/optrove/sif/raw/HEAD/ARWHEAD.SIF)
# see also (https://vanderbei.princeton.edu/ampl/nlmodels/cute/arwhead.mod)
#  classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function OARWHEADFun(x::Vector{<:Real})
    # sum {i in 1..N-1} (-4*x[i]+3.0) + sum {i in 1..N-1} (x[i]^2+x[N]^2)^2
    n = length(x)
    sum((-4 * x[i] + 3.0) for i in 1:n-1) + sum((x[i]^2 + x[n]^2)^2 for i in 1:n-1)
end
function OARWHEADGrad(x::Vector{<:Real})
    vcat(-4 .+ 4 * x[1:end-1] .* (x[1:end-1] .^ 2 .+ x[end]^2),
        4 * x[end] * sum(x[1:end-1] .^ 2 .+ x[end]^2))
end


# name: PENALTY1 see (https://bitbucket.org/optrove/sif/raw/HEAD/PENALTY1.SIF) 
# see also (https://vanderbei.princeton.edu/ampl/nlmodels/cute/penalty1.mod)
# classification SUR2-AN-V-0
# S: the objective function Ois a sum of squares
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: the problem is academic, that is, has been constructed specifically by researchers to test one or more algorithms,
# N: the problem description does not contain any explicit internal variables.
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.

function OPENALTY1Fun(x::Vector{<:Real})
    a = 1e-5
    N = length(x)
    # sum {i in 1..N} a*(x[i]-1)^2 + ( sum {j in 1..N} x[j]^2 - 1/4 )^2;
    sum(a * (x[i] - 1)^2 for i in 1:N) + (sum(x[i]^2 for i in 1:N) - 0.25)^2
end
function OPENALTY1Grad(x::Vector{<:Real})
    t = sum(x .^ 2)
    c = 1e-5
    2 * c .* (x .- 1) + 4 * (t - 0.25) .* x
end

# NAME: DIXON3DQ see (https://bitbucket.org/optrove/sif/raw/HEAD/DIXON3DQ.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/dixon3dq.mod)
# classification QUR2-AN-V-0
# Q: the objective function Ois quadratic,
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function ODIXON3DQFun(x::Vector{<:Real})
    # (x[1]-1.0)^2 + sum {j in 2..n-1} (x[j]-x[j+1])^2 + (x[n]-1.0)^2
    n = length(x)
    (x[1] - 1.0)^2 + sum((x[j] - x[j+1])^2 for j in 2:n-1) + (x[n] - 1.0)^2
end
function ODIXON3DQGrad(x::Vector{<:Real})
    n = length(x)
    return vcat(2 * (x[1] - 1),
        2 * (x[2] - x[3]),
        2 * (2 * x[3:n-1] - x[2:n-2] - x[4:n]),
        2 * (2 * x[n] - x[n-1] - 1)
    )
end

# NAME: GENHUMPS, see (https://bitbucket.org/optrove/sif/raw/HEAD/GENHUMPS.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/genhumps.mod)
# classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function OGENHUMPSFun(x::Vector{<:Real}; ζ=20)
    # sum {i in 1..N-1} ( sin (zeta*x[i])^2*sin(zeta*x[i+1])^2+0.05*(x[i]^2+x[i+1]^2) );
    N = length(x)
    sum((sin(ζ * x[i]) * sin(ζ * x[i+1]))^2 + 0.05 * (x[i]^2 + x[i+1]^2) for i in 1:(N-1))
end
function OGENHUMPSGrad(x::Vector{<:Real}; ζ=20)
    TNTHX = 0.1 * x
    ZX = ζ * x
    ZX2 = 2 * ZX
    SZX = (sin.(ZX)) .^ 2
    SZX2 = sin.(ZX2)
    ZSZX2 = ζ * SZX2
    vcat(
        ZSZX2[1] * SZX[2] + TNTHX[1],
        ZSZX2[2:end-1] .* (SZX[1:end-2] .+ SZX[3:end]) .+ (2 * TNTHX[2:end-1]),
        ZSZX2[end] * SZX[end-1] + TNTHX[end],
    )
end

# NAME: ENGVAL1 see (https://bitbucket.org/optrove/sif/raw/HEAD/ENGVAL1.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/engval1.mod)
# classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function OENGVAL1Fun(x::Vector{<:Real})
    # sum {i in 1..N-1} (x[i]^2+x[i+1]^2)^2 + sum {i in 1..N-1} (-4*x[i]+3.0);
    N = length(x)
    sum((x[i]^2 + x[i+1]^2)^2 for i in 1:(N-1)) + sum(-4 * x[i] + 3.0 for i in 1:(N-1))
end
function OENGVAL1Grad(x::Vector{<:Real})
    X2 = x .^ 2
    vcat(
        4 * (x[1] * (X2[1] + X2[2]) - 1),
        4 * (x[2:end-1] .* (X2[1:end-2] + 2 * X2[2:end-1] + X2[3:end]) .- 1),
        4 * x[end] * (X2[end-1] + X2[end])
    )
end
# NAME: DIXMAANH see (https://bitbucket.org/optrove/sif/raw/HEAD/DIXMAANH.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/dixmaanh.mod)
# classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function ODIXMAANHFun(x::Vector{<:Real}; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
    N = length(x)
    if N % 3 != 0
        throw("The length of `x0` must be divisible by 3")
    end
    M = fld(N, 3)
    1.0 + α * sum(x[i]^2 * (i / N)^K[1] for i in 1:N) +
    β * sum(x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i in 1:N-1) +
    γ * sum(x[i]^2 * x[i+M]^4 for i in 1:2*M) +
    δ * sum(x[i] * x[i+2*M] * (i / N)^K[4] for i in 1:M)
end
function ODIXMAANHGrad(x::Vector{<:Real}; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
    N = length(x)
    if N % 3 != 0
        throw("The length of `x0` must be divisible by 3")
    end
    M = fld(N, 3)
    X2 = x .^ 2
    X4 = X2 .* X2
    IOverN = [(i / N)^K[1] for i in 1:N]
    IOverM = [(i / N)^K[4] for i in 1:M]
    TwoAlphaXOverN = 2 * α * IOverN .* x
    OnePlusX = 1 .+ x
    OnePlus2X = OnePlusX .+ x
    OnePlusX2 = OnePlusX .^ 2
    G1 = TwoAlphaXOverN[1] + 2 * β * x[1] * X2[2] * OnePlusX2[2] + 2 * γ * x[1] * X4[1+M] + δ * IOverM[1] * x[1+2*M]
    G2M = TwoAlphaXOverN[2:M] .+ 2 * β * (X2[1:M-1] .* x[2:M] .* OnePlusX[2:M] .* OnePlus2X[2:M] .+ x[2:M] .* X2[3:M+1] .* OnePlusX2[3:M+1]) .+
          2 * γ * x[2:M] .* X4[2+M:2*M] .+ δ * IOverM[2:M] .* x[2+2*M:3*M]

    # check
    GM12M = TwoAlphaXOverN[M+1:2*M] .+
            2 * β * (
                X2[M:2*M-1] .* x[M+1:2*M] .* OnePlusX[M+1:2*M] .* OnePlus2X[M+1:2*M] .+ x[M+1:2*M] .* X2[M+2:2*M+1] .* OnePlusX2[M+2:2*M+1]
            ) .+
            γ * (2 * x[M+1:2*M] .* X4[2*M+1:3*M] .+ 4 * X2[1:M] .* x[M+1:2*M] .* X2[M+1:2*M])


    # check
    G2M1Nm1 = TwoAlphaXOverN[2*M+1:N-1] .+
              2 * β * (
                  X2[2*M:N-2] .* x[2*M+1:N-1] .* OnePlusX[2*M+1:N-1] .* OnePlus2X[2*M+1:N-1] .+ x[2*M+1:N-1] .* X2[2*M+2:N] .* OnePlusX2[2*M+2:N]
              ) .+
              (γ * 4 * X2[M+1:2*M-1] .* x[2*M+1:N-1] .* X2[2*M+1:N-1]) .+
              δ * IOverM[1:M-1] .* x[1:M-1]

    GN = TwoAlphaXOverN[N] + 2 * β * (X2[N-1] * x[N] * OnePlusX[N] * OnePlus2X[N]) +
         γ * (4 * X2[2*M] * x[N] * X2[N]) + δ * IOverM[M] .* x[M]
    vcat(
        G1,
        G2M,
        GM12M,
        G2M1Nm1,
        GN
    )

end

function ODIXMAANIFun(x::Vector{<:Real})
    DIXMAANHFun(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
end
function ODIXMAANIGrad(x::Vector{<:Real})
    DIXMAANHGrad(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
end





## new from old with exact projections
# 
# 19 Χ*****************************************Χ
function COExponetialIII(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(0.0, Inf64)), x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)), n)
    NECProblem(x -> OExponetialIII(x), Ω)
end

"""
Problem 17
"""
# 20
function COPolynomialI(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(0.0, Inf64)), x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)), n)
    NECProblem(x -> OPolynomialI(x), Ω)
end


"""
Problem 18
"""
# 21
function COSmoothSine(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(0.0, Inf64)), x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)), n)
    NECProblem(x -> OSmoothSine.(x), Ω)
end


"""
Problem 19
"""
# 22
function COPolynomialSineCosine(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(0.0, Inf64)), x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)), n)
    NECProblem(x -> OPolynomialSineCosine(x), Ω)
end


"""
Problem 20
"""
# 23
function COExponetialI(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(0.0, Inf64)), x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)), n)
    NECProblem(x -> OExponetialI.(x), Ω)
end


"""
Problem 21
"""
# 24
function CONonsmoothSine(n::Int)
    Ω = Polyhedral(x -> projectOnTriangle(x; bounds=(0, Inf64), β=n), x -> projectOnTriangleCheck(x; bounds=(0, Inf64), β=n), n)
    NECProblem(x -> ONonsmoothSine.(x), Ω)
end


"""
Problem 22
"""
# 25
function COModifiedNonsmoothSine(n::Int)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=(-1, Inf64)), x -> projectOnTriangleCheck(x; β=n, bounds=(-1, Inf64)), n)
    NECProblem(x -> OModifiedNonsmoothSine.(x), Ω)
end

"""
Problem 23
"""
# 26
function COModifiedNonsmoothSine2(n::Int)
    Ω = Polyhedral(x -> projectOnTriangle(x; bounds=(-1, Inf64), β=n), x -> projectOnTriangleCheck(x; bounds=(-1, Inf64), β=n), n)
    NECProblem(x -> OModifiedNonsmoothSine2.(x), Ω)
end


"""
Problem 24
"""
# 27
function COExponetialSineCosine(n::Int)
    Ω = Polyhedral(x -> projectOnTriangle(x; bounds=(-1, Inf64), β=n), x -> projectOnTriangleCheck(x; bounds=(-1, Inf64), β=n), n)
    NECProblem(x -> OExponetialSineCosine.(x), Ω)
end

"""
Problem 25
"""
# 28
function COModifiedTrigI(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(-3.0, Inf64)), x -> projectOnBoxCheck(x; bounds=(-3.0, Inf64)), n)
    NECProblem(x -> OModifiedTrigI(x), Ω)
end


# 29 (not used)
function COModifiedTridiagonal(n::Int)
    Ω = Polyhedral(x -> x, x -> true, n)
    NECProblem(x -> OModifiedTridiagonalV2(x), Ω)
end


"""
Problem 26
"""
# 30 (29)
function COLogarithmic(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(-1.0, Inf64)), x -> projectOnBoxCheck(x; bounds=(-1.0, Inf64)), n)
    NECProblem(x -> OLogarithmic(x), Ω)
end


"""
Problem 27
"""
# 31 (30)
function CONonmoothLogarithmic(n::Int)
    Ω = Polyhedral(x -> projectOnBox(x; bounds=(0, Inf64)), x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)), n)
    NECProblem(x -> ONonmoothLogarithmic(x), Ω)
end

# 32 (31) not used
function COARWHEADGrad(n::Int)
    Ω = Polyhedral(x -> projectOnTriangle(x; bounds=(0, Inf64), β=n), x -> projectOnTriangleCheck(x; bounds=(0, Inf64), β=n), n)
    NECProblem(x -> OARWHEADGrad(x), Ω)
end


"""
Problem 28
"""
# 33 (32)
function COENGVAL1Grad(n::Int)
    Ω = Polyhedral(x -> x, x -> true, n)
    NECProblem(x -> OENGVAL1Grad(x), Ω)
end

