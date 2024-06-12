"""
Source:
Problem 4.6 in
Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, 
Modified optimal Perry conjugate gradient method for solving system of monotone equations with applications,Applied Numerical Mathematics,
Volume 184, 2023, Pages 431-445, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2022.10.016.

Projected on [0, ∞] or [-1, ∞]

"""
function PolynomialSineCosine(x)
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
function ExponetialI(x)
    exp(x) - 1
end



"""
Source:
Problem 3 in Appendix
William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

"""
function ExponetialIII(x)
    X2 = x .^ 2
    Ii = [i / 10 for i in 1:length(x)]
    Y = vcat(
        X2[1:end-1] .+ exp.(-X2[1:end-1]),
        exp(-X2[end])
    )
    Ii .* (1 .- Y)
end


"""
Source:
Problem 9 
J. Sabi’u, A. Shah, M.Y. Waziri, M.K. Dauda, A new hybrid approach for solving large-scale monotone nonlinear equations, 
J. Math. Fund. Sci. 52 (2020) 17–26 https://doi.org/10.5614/j.math.fund.sci.2020.52.1.2


"""
function PolynomialI(x)
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
W.J. Zhou, D.H. Li, A globally convergent BFGS method for nonlinear monotone equations without any merit functions, 
Math. Comput. 77 (264) (2008) 2231–2240.
Projected on [-2, ∞]


"""
function SmoothSine(x)
    v = try

        2 * x - sin(x)
    catch
        println(max(x...), min(x...))
        throw("problem with this z")
    end
    v
end



"""
Source:
Problem 1 in
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.tional and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

Project to Rn₊
"""
function NonsmoothSine(x)
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
function ModifiedNonsmoothSine(x)
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
function ModifiedNonsmoothSine2(x)
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
function ExponetialSineCosine(x)
    (exp.(x)) .^ 2 + 3 * sin.(x) * cos.(x) - 1
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
function ModifiedTrigI(x)
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
function Tridiagonal(x)
    # Ii = [i for i in 1:length(x)]
    x0 = vcat(0, x)
    x1 = vcat(x, 0)
    h = 1 / (1 + length(x))
    x .- exp.(cos.(h .* (x0[1:end-1] .+ x .+ x1[2:end])))
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
function ModifiedTridiagonal(x)
    n = length(x)
    vcat(
        x[1] - exp(cos((x[1] + x[2]) / (n + 1))),
        -x[2:end-1] - exp.(cos.((x[1:end-2] .+ x[2:end-1] .+ x[3:end]) ./ (n + 1))),
        x[end] - exp(cos((x[end-1] + x[end]) / (n + 1)))
    )
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
function Logarithmic(x)
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
function NonmoothLogarithmic(x)
    log.(abs.(x) .+ 1) .- (x / length(x))
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
function ARWHEADFun(x)
    # sum {i in 1..N-1} (-4*x[i]+3.0) + sum {i in 1..N-1} (x[i]^2+x[N]^2)^2
    n = length(x)
    sum((-4 * x[i] + 3.0) for i in 1:n-1) + sum((x[i]^2 + x[n]^2)^2 for i in 1:n-1)
end
function ARWHEADGrad(x)
    vcat(-4 .+ 4 * x[1:end-1] .* (x[1:end-1] .^ 2 .+ x[end]^2),
        4 * x[end] * sum(x[1:end-1] .^ 2 .+ x[end]^2))
end


# name: PENALTY1 see (https://bitbucket.org/optrove/sif/raw/HEAD/PENALTY1.SIF) 
# see also (https://vanderbei.princeton.edu/ampl/nlmodels/cute/penalty1.mod)
# classification SUR2-AN-V-0
# S: the objective function is a sum of squares
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: the problem is academic, that is, has been constructed specifically by researchers to test one or more algorithms,
# N: the problem description does not contain any explicit internal variables.
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.

function PENALTY1Fun(x)
    a = 1e-5
    N = length(x)
    # sum {i in 1..N} a*(x[i]-1)^2 + ( sum {j in 1..N} x[j]^2 - 1/4 )^2;
    sum(a * (x[i] - 1)^2 for i in 1:N) + (sum(x[i]^2 for i in 1:N) - 0.25)^2
end
function PENALTY1Grad(x)
    t = sum(x .^ 2)
    c = 1e-5
    2 * c .* (x .- 1) + 4 * (t - 0.25) .* x
end

# NAME: DIXON3DQ see (https://bitbucket.org/optrove/sif/raw/HEAD/DIXON3DQ.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/dixon3dq.mod)
# classification QUR2-AN-V-0
# Q: the objective function is quadratic,
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function DIXON3DQFun(x)
    # (x[1]-1.0)^2 + sum {j in 2..n-1} (x[j]-x[j+1])^2 + (x[n]-1.0)^2
    n = length(x)
    (x[1] - 1.0)^2 + sum((x[j] - x[j+1])^2 for j in 2:n-1) + (x[n] - 1.0)^2
end
function DIXON3DQGrad(x)
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
function GENHUMPSFun(x; ζ=20)
    # sum {i in 1..N-1} ( sin (zeta*x[i])^2*sin(zeta*x[i+1])^2+0.05*(x[i]^2+x[i+1]^2) );
    N = length(x)
    sum((sin(ζ * x[i]) * sin(ζ * x[i+1]))^2 + 0.05 * (x[i]^2 + x[i+1]^2) for i in 1:(N-1))
end
function GENHUMPSGrad(x; ζ=20)
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
function ENGVAL1Fun(x)
    # sum {i in 1..N-1} (x[i]^2+x[i+1]^2)^2 + sum {i in 1..N-1} (-4*x[i]+3.0);
    N = length(x)
    sum((x[i]^2 + x[i+1]^2)^2 for i in 1:(N-1)) + sum(-4 * x[i] + 3.0 for i in 1:(N-1))
end
function ENGVAL1Grad(x)
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
function DIXMAANHFun(x; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
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
function DIXMAANHGrad(x; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
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

function DIXMAANIFun(x)
    DIXMAANHFun(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
end
function DIXMAANIGrad(x)
    DIXMAANHGrad(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
end
"""
Source:
Problem 2 in Appendix

William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

"""
function ExponetialII(x)
    a = [i / 10 for i in 2:length(x)]
    vcat(exp(x[1] - 1), a .* exp.(x[2:end]) .+ x[2:end] .- 1)
end

"""
Source:
Problem 17 in Appendix (modified)

William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

"""
function ExponetialIV(x)
    n = length(x)
    a = collect(1:length(x))
    (a .* (exp.(x)) / n) .- 1
end


"""
Source:
Problem 2 
Wang, C., Wang, Y. & Xu, C. 
A projection method for a system of nonlinear monotone equations with convex 
constraints. Math Meth Oper Res 66, 33–46 (2007). 
https://doi.org/10.1007/s00186-006-0140-y

"""
function ArcTan(x)
    n = length(x)
    d = Normal(0.0, Float64(n))
    M = rand(d, n, n)
    # M = 5 * [0.33435 0.535161 0.561637 0.344291 0.111954
    #     0.360864 0.508599 0.821533 0.694135 0.953729
    #     0.876373 0.747615 0.336939 0.0355163 0.558138
    #     0.983291 0.355772 0.209947 0.251682 0.410537
    #     0.270499 0.514839 0.198478 0.883239 0.74523]
    M = 0.5 * (M - M')
    q = rand(d, n)
    # q = [0.7845171302827574
    #     0.5478931492840216
    #     0.8419546163518192
    #     0.3928429699008852
    #     0.21783349212498404]
    ρ = 100
    q0 = 1.5 * ones(n)
    F(y) = ρ * atan.(y .- 2) + M * y + q
    F(x) - F(q0)
end
"""
Source:
Problem 3 
Wang, C., Wang, Y. & Xu, C. 
A projection method for a system of nonlinear monotone equations with convex 
constraints. Math Meth Oper Res 66, 33–46 (2007). 
https://doi.org/10.1007/s00186-006-0140-y

"""
function ArcTan2(x)
    n = length(x)
    d1 = Normal(0.0, 0.5)
    d2 = Normal(0.0, 250)
    A = rand(d1, n, n)
    B = rand(d1, n, n)
    B = 0.5 * (B - B')
    M = A * A' + B
    q = rand(d2, n)
    q02 = 0.9 * ones(n)
    a = rand(1:100, n)
    D(y) = a .* atan.(y)
    F(y) = D(y) + M * y + q
    F(x) - F(q02)
end

"""
Source:
Problem 4 
Wang, C., Wang, Y. & Xu, C. 
A projection method for a system of nonlinear monotone equations with convex 
constraints. Math Meth Oper Res 66, 33–46 (2007). 
https://doi.org/10.1007/s00186-006-0140-y
"""
function StrictlyMonotone(x)
    A = [
        1 0 0 0
        0 1 -1 0
        0 1 1 0
        0 0 0 0
    ]
    b = [-10; 1; -3; 0]
    A * x + [(x[1])^3; (x[2])^3; 2(x[3])^3; 2(x[4])^3] + b
end

"""
Source:
Problem 1 in 

Abubakar, A. B., Kumam, P., & Mohammad, H. (2020). 
A note on the spectral gradient projection method for nonlinear monotone equations with applications. 
Computational and Applied Mathematics, 39(2), 129. https://doi.org/10.1007/s40314-020-01151-5

"""
function ExponetialIIModified(x)
    vcat(exp(x[1]) - 1, exp.(x[2:end]) .+ x[2:end] .- 1)
end

"""
Source:
Problem 4.2 in 

Liu, J., & Feng, Y. (2019). 
A derivative-free iterative method for nonlinear monotone equations with convex constraints. 
Numerical Algorithms, 82(1), 245–262. https://doi.org/10.1007/s11075-018-0603-2

"""
function DiscreteBV(x::Vector{Float64})
    n = length(x)
    h = 1 / (n + 1)
    F = zeros(n)

    # Equation 1
    F[1] = 2 * x[1] + 0.5 * h^2 * (x[1] + h)^3 - x[2]

    # Equations 2 to n-1
    for i in 2:n-1
        F[i] = 2 * x[i] + 0.5 * h^2 * (x[i] + i * h)^3 - x[i-1] + x[i+1]
    end

    # Equation n
    F[n] = 2 * x[n] + 0.5 * h^2 * (x[n] + n * h)^3 - x[n-1]

    return F
end

"""
Source:
Problem 4.3 in 

Liu, J., & Feng, Y. (2019). 
A derivative-free iterative method for nonlinear monotone equations with convex constraints. 
Numerical Algorithms, 82(1), 245–262. https://doi.org/10.1007/s11075-018-0603-2

"""
function TriExponential(x::Vector{Float64})
    n = length(x)
    F = zeros(Float64, n)

    # Equation 1
    F[1] = 3 * x[1]^3 + 2 * x[2] - 5 + sin(x[1] - x[2]) * sin(x[1] + x[2])

    # Equations 2 to n-1
    for i in 2:n-1
        F[i] = -x[i-1] * exp(x[i-1] - x[i]) + x[i] * (4 + 3 * x[i]^2) +
               2 * x[i+1] + sin(x[i] - x[i+1]) * sin(x[i] + x[i-1]) - 8
    end

    # Equation n
    F[n] = -x[n-1] * exp(x[n-1] - x[n]) + 4 * x[n] - 3

    return F
end


"""
Source:
Problem 4.7 in 

Liu, J., & Feng, Y. (2019). 
A derivative-free iterative method for nonlinear monotone equations with convex constraints. 
Numerical Algorithms, 82(1), 245–262. https://doi.org/10.1007/s11075-018-0603-2

"""
function PolynomialModified(x::Vector{Float64})
    n = length(x)
    F = zeros(Float64, n)

    # Equation 1
    F[1] = 2.5 * x[1] + x[2] - 1

    # Equations 2 to n-1
    for i in 2:n-1
        F[i] = x[i-1] + 2.5 * x[i] + x[i+1] - 1
    end

    # Equation n
    F[n] = x[n-1] + 2.5 * x[n] - 1

    return F
end

"""
Source:
Problem 1 in 
Cruz, W. L. (2017). 
A spectral algorithm for large-scale systems of nonlinear monotone equations. 
Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

"""
function MinimizeMaximizeFunction(x::Vector{Float64})
    n = length(x)
    F = zeros(Float64, n)

    for i in 1:n
        F[i] = min(min(abs(x[i]), x[i]^2), max(abs(x[i]), x[i]^3))
    end

    return F
end

"""
Source:
Problem 4 in 
Cruz, W. L. (2017). 
A spectral algorithm for large-scale systems of nonlinear monotone equations. 
Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

"""
function SinusoidalFunction(x::Vector{Float64})
    n = length(x)
    F = zeros(Float64, n)

    F[1] = 2 * x[1] + sin(x[1]) - 1

    for i in 2:n-1
        F[i] = -2 * x[i-1] + 2 * x[i] + sin(x[i]) - 1
    end

    F[n] = 2 * x[n] + sin(x[n]) - 1

    return F
end



"""
Source:
Problem 7 in 
Cruz, W. L. (2017). 
A spectral algorithm for large-scale systems of nonlinear monotone equations. 
Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

"""
function LinearSystemMatrix(x::Vector{Float64})
    n = length(x)
    A = (5 / 2) * I(n) + diagm(-1 => ones(n - 1), 1 => ones(n - 1))
    b = -ones(Float64, n)

    return A * x + b
end

"""
Source:
Problem 8 in 
Cruz, W. L. (2017). 
A spectral algorithm for large-scale systems of nonlinear monotone equations. 
Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

"""
function CubicPolynomial(x::Vector{Float64})
    n = length(x)
    F = zeros(Float64, n)

    F[1] = (1 / 3) * x[1]^3 + (1 / 2) * x[2]^2

    for i in 2:n-1
        F[i] = -(1 / 2) * x[i]^2 + (i / 3) * x[i]^3 + (1 / 2) * x[i+1]^2
    end

    F[n] = -(1 / 2) * x[n]^2 + (n / 3) * x[n]^3

    return F
end

"""
Source:
Problem 9 in 
Cruz, W. L. (2017). 
A spectral algorithm for large-scale systems of nonlinear monotone equations. 
Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

"""
function QuadraticSumFunction(x::Vector{Float64})
    n = length(x)
    F = zeros(Float64, n)
    smx = sum(x)
    for i in 1:n
        F[i] = x[i] + ((smx - x[i]^2) / n) + i
    end

    return F
end

# ==================================================================================================================== #
"""
=== Problems in Gonçalves, M. L. N., & Menezes, T. C. (2023) ===
"""

"""
Problem 1 in

Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

See ExponetialI

"""
function ConstrainedExponetialI(n::Int)
    bounds = (-1, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=n, bvalue=n)
    NECProblem(x -> ExponetialI.(x), Ω)
end

"""
Problem 2 in

Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

See ModifiedNonsmoothSine
"""
function ConstrainedModifiedNonsmoothSine(n::Int)
    # β is bvalue
    bounds = (-1, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=n, bvalue=n)
    NECProblem(x -> ModifiedNonsmoothSine.(x), Ω)
end

"""
Source:
Problem 3 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See ArcTan
"""
function ConstrainedArcTan(n::Int)
    bounds = (0, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n; lower=0, bvalue=2n, upper=n)
    return NECProblem(x -> ArcTan(x), Ω)
end

"""
Source:

Problem 4 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See ArcTan2
"""
function ConstrainedArcTan2(n::Int)
    bounds = (0, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; bounds=bounds),
        x -> projectOnTriangleCheck(x; bounds=bounds), n)
    # Ω = createPolyhedral(n; lower=0, bvalue=n, upper=n)
    return NECProblem(x -> ArcTan2(x), Ω)
end


"""
Problem 5 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See StrictlyMonotone
"""
function ConstrainedStrictlyMonotone(n::Int)
    bounds = (-1, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=4, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=4, bounds=bounds), n)
    # Ω = createPolyhedral(n; lower=-1, bvalue=4, upper=n)
    return NECProblem(x -> StrictlyMonotone(x), Ω)
end

"""
Source:
Problem 6 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z
    
    See ExponetialII
"""
function ConstrainedExponetialII(n::Int)
    bounds = (-1, 2)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=2, bvalue=n)
    return NECProblem(x -> ExponetialIIModified(x), Ω)
end

"""
Source:
Problem 7 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See Logarithmic
"""
function ConstrainedLogarithmic(n::Int)
    bounds = (-1, 2)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=2, bvalue=n)
    return NECProblem(x -> Logarithmic(x), Ω)
end


"""
Source:
Problem 8 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See NonsmoothSine    
"""
function ConstrainedNonsmoothSine(n::Int)
    bounds = (0, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=0, upper=n, bvalue=n)
    return NECProblem(x -> NonsmoothSine.(x), Ω)
end


"""
Source:
Problem 9 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See ExponetialIV
"""

function ConstrainedExponetialIV(n::Int)
    bounds = (-1, 7)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=1.1 * n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=1.1 * n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=7, bvalue=1.1 * n)
    return NECProblem(x -> ExponetialIV(x), Ω)
end

"""
Source:
Problem 10 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See Tridiagonal
"""

function ConstrainedTridiagonal(n::Int)
    bounds = (0, exp(1))
    Ω = Polyhedral(x -> projectOnTriangle(x; β=exp(1) * n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=exp(1) * n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=0, upper=exp(1), bvalue=exp(1) * n)
    return NECProblem(x -> Tridiagonal(x), Ω)
end

"""
Χ*****************************************Χ
Source:
Problem 11 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z
    
"""
function ConstrainedMinimizeMaximizeFunction(n::Int)
    Ω = createPolyhedral(n, 10; lower=-1, upper=2, bvalue=n)
    return NECProblem(x -> MinimizeMaximizeFunction(x), Ω)
end

"""
Χ*****************************************Χ
Source:
Problem 12 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z
    
"""
function ConstrainedSinusoidalFunction(n::Int)
    bounds = (-1, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=20 * n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=20 * n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=n, bvalue=20 * n)
    return NECProblem(x -> SinusoidalFunction(x), Ω)
end

"""
Problem 11
Source:
Problem 13 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z
    
"""
function ConstrainedLinearSystemMatrix(n::Int)
    bounds = (-1, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=n, bvalue=n)
    return NECProblem(x -> LinearSystemMatrix(x), Ω)
end


"""
Problem 12
Source:
Problem 14 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z
    
"""
function ConstrainedCubicPolynomial(n::Int)
    bounds = (-n, 1)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-n, upper=1, bvalue=n)
    return NECProblem(x -> CubicPolynomial(x), Ω)
end


"""
Problem 13
Source:
Problem 15 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z
    
"""
function ConstrainedQuadraticSumFunction(n::Int)
    bounds = (-n, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-n, upper=n, bvalue=n)
    return NECProblem(x -> QuadraticSumFunction(x), Ω)
end

"""
Problem 14
Source:
Problem 16 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See Tridiagonal
"""

function ConstrainedDiscreteBV(n::Int)
    bounds = (-n, 1)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=1, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=1, bounds=bounds), n)

    # Ω = createPolyhedral(n, lower=-n, upper=1, bvalue=1)
    return NECProblem(x -> DiscreteBV(x), Ω)
end

"""
Problem 15
Source:
Problem 17 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See Tridiagonal
"""

function ConstrainedTriExponential(n::Int)
    bounds = (-1, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=n, bvalue=n)
    return NECProblem(x -> TriExponential(x), Ω)
end

"""
Problem 16
Source:
Problem 18 in 
Gonçalves, M. L. N., & Menezes, T. C. (2023). 
A framework for convex-constrained monotone nonlinear equations and its special cases. Computational and Applied Mathematics, 42(7), 306. https://doi.org/10.1007/s40314-023-02446-z

    See Tridiagonal
"""

function ConstrainedPolynomialModified(n::Int)
    bounds = (-1, n)
    Ω = Polyhedral(x -> projectOnTriangle(x; β=n, bounds=bounds),
        x -> projectOnTriangleCheck(x; β=n, bounds=bounds), n)
    # Ω = createPolyhedral(n, lower=-1, upper=n, bvalue=n)
    return NECProblem(x -> PolynomialModified(x), Ω)
end

