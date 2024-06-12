function init_experiment(n::Int)
    initial_pts = createStartingPoints()
    # :tenth => n -> 0.1 * ones(n),
    # :negativetenth => n -> -0.1 * ones(n),
    # :allones => n -> ones(n),
    # :negativeones => n -> -ones(n),
    # :ntenthn => n -> begin
    #     i = (n - 1) / n
    #     vcat(i, 0.1 * ones(n - 2), i)
    # end,
    # :allzeros => n -> zeros(n),
    # :tenzeros => n -> vcat(10, zeros(n - 1)),
    # :ninezeros => n -> vcat(9, zeros(n - 1)),
    # :threezeros => n -> vcat(3, zeros(n - 1)),
    # :threeszeros => n -> [iseven(i) ? 0 : 3 for i in 1:n],
    # :zerotwos => n -> vcat(0, 2 * ones(n - 1)),
    # :zeroones => n -> vcat(0, ones(n - 1)),
    # :oneszero => n -> vcat(ones(n - 1), 0),
    # :halfn => n -> [1 / 2^i for i in 1:n],
    # :negativehalfn => n -> -[1 / 2^i for i in 1:n],
    # :nth => n -> [1 / i for i in 1:n],
    # :negativenth => n -> -[1 / i for i in 1:n]
    p1 = ConstrainedExponetialI(n)
    p1_x0_names = vcat(:tenth, :allones, :ntenthn, :negativeones)
    x0_p1 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p1_x0_names)

    p2 = ConstrainedModifiedNonsmoothSine(n)
    p2_x0_names = vcat(:tenth, :allones, :allzeros, :negativeones)
    x0_p2 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p2_x0_names)

    p3 = ConstrainedArcTan(5)
    p3_x0_names = vcat(:tenzeros, :ninezeros, :threeszeros, :zerotwos)
    x0_p3 = map(x -> (x, get(initial_pts, x, m -> ones(m))(5)), p3_x0_names)

    p4 = ConstrainedArcTan2(10)
    p4_x0_names = vcat(:allones, :tenth, :halfn, :nth)
    x0_p4 = map(x -> (x, get(initial_pts, x, m -> ones(m))(10)), p4_x0_names)
    p_like_4_x0_names = vcat(:allones, :tenth, :ntenthn, :nth)
    x0_p_like_4 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p_like_4_x0_names)

    p5 = ConstrainedStrictlyMonotone(4)
    p5_x0_names = vcat(:allzeros, :threezeros, :oneszero, :zeroones)
    x0_p5 = map(x -> (x, get(initial_pts, x, m -> ones(m))(4)), p5_x0_names)

    p6 = ConstrainedExponetialII(n)
    p6_x0_names = copy(p_like_4_x0_names)
    x0_p6 = copy(x0_p_like_4)

    p7 = ConstrainedLogarithmic(n)
    p7_x0_names = vcat(:allones, :tenth, :halfn, :nth)
    x0_p7 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p7_x0_names)

    p8 = ConstrainedNonsmoothSine(n)
    p8_x0_names = copy(p_like_4_x0_names)
    x0_p8 = copy(x0_p_like_4)

    p9 = ConstrainedExponetialIV(n)
    p9_x0_names = copy(p_like_4_x0_names)
    x0_p9 = copy(x0_p_like_4)

    p10 = ConstrainedTridiagonal(n)
    p10_x0_names = copy(p_like_4_x0_names)
    x0_p10 = copy(x0_p_like_4)

    p11 = ConstrainedMinimizeMaximizeFunction(n)
    p11_x0_names = vcat(:allones, :onestimes8, :negativeones, :onestimes2)
    x0_p11 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p11_x0_names)

    p12 = ConstrainedSinusoidalFunction(n)
    p12_x0_names = copy(p_like_4_x0_names)
    x0_p12 = copy(x0_p_like_4)

    p13 = ConstrainedLinearSystemMatrix(n)
    p13_x0_names = copy(p_like_4_x0_names)
    x0_p13 = copy(x0_p_like_4)

    p14 = ConstrainedCubicPolynomial(n)
    p14_x0_names = vcat(:negativeones, :negativetenth, :negativehalfn, :negativenth)
    x0_p14 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p14_x0_names)

    p15 = ConstrainedQuadraticSumFunction(n)
    p15_x0_names = copy(p_like_4_x0_names)
    x0_p15 = copy(x0_p_like_4)

    p16 = ConstrainedDiscreteBV(n)
    p16_x0_names = copy(p4_x0_names)
    x0_p16 = copy(x0_p14)

    p17 = ConstrainedTriExponential(n)
    # x0_p17 = copy(x0_p_like_4)
    p17_x0_names = vcat(:ntenthn, :tenth, :halfn, :nth)
    x0_p17 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p17_x0_names)

    p18 = ConstrainedPolynomialModified(n)
    p18_x0_names = copy(p_like_4_x0_names)
    x0_p18 = copy(x0_p_like_4)


    ## old problems with
    p19 = COExponetialIII(n)
    p19_x0_names = copy(p_like_4_x0_names)
    x0_p19 = copy(x0_p_like_4)

    p20 = COPolynomialI(n)
    p20_x0_names = copy(p_like_4_x0_names)
    x0_p20 = copy(x0_p_like_4)

    p21 = COSmoothSine(n)
    p21_x0_names = copy(p_like_4_x0_names)
    x0_p21 = copy(x0_p_like_4)


    p22 = COPolynomialSineCosine(n)
    p22_x0_names = vcat(:allones, :onestimes8, :threeszeros, :halfn)
    x0_p22 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p22_x0_names)

    p23 = COExponetialI(n)
    p23_x0_names = copy(p_like_4_x0_names)
    x0_p23 = copy(x0_p_like_4)

    p24 = CONonsmoothSine(n)
    p24_x0_names = copy(p_like_4_x0_names)
    x0_p24 = copy(x0_p_like_4)


    p25 = COModifiedNonsmoothSine(n)
    p25_x0_names = copy(p_like_4_x0_names)
    x0_p25 = copy(x0_p_like_4)


    p26 = COModifiedNonsmoothSine2(n)
    p26_x0_names = copy(p_like_4_x0_names)
    x0_p26 = copy(x0_p_like_4)

    p27 = COExponetialSineCosine(n)
    p27_x0_names = copy(p_like_4_x0_names)
    x0_p27 = copy(x0_p_like_4)

    p28 = COModifiedTrigI(n)
    p28_x0_names = copy(p_like_4_x0_names)
    x0_p28 = copy(x0_p_like_4)


    p29 = COLogarithmic(n)
    p29_x0_names = vcat(:allones, :tenth, :halfn, :nth)
    x0_p29 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), p29_x0_names)


    p30 = CONonmoothLogarithmic(n)
    p30_x0_names = copy(p_like_4_x0_names)
    x0_p30 = copy(x0_p_like_4)

    p31 = COARWHEADGrad(n)
    p31_x0_names = copy(p_like_4_x0_names)
    x0_p31 = copy(x0_p_like_4)

    p32 = COENGVAL1Grad(n)
    p32_x0_names = copy(p_like_4_x0_names)
    x0_p32 = copy(x0_p_like_4)

    # x0_p33 = map(x -> (x, get(initial_pts, x, m -> ones(m))(n)), vcat(:allones, :onestimeshalf, :onestimes8, :onestimes11half))
    # p33 = COModifiedTridiagonal(n)

    [
        ("p1", p1, x0_p1, p1_x0_names),
        ("p2", p2, x0_p2, p2_x0_names),
        ("p3", p3, x0_p3, p3_x0_names),
        ("p4", p4, x0_p4, p4_x0_names),
        ("p5", p5, x0_p5, p5_x0_names),
        ("p6", p6, x0_p6, p6_x0_names),
        ("p7", p7, x0_p7, p7_x0_names),
        ("p8", p8, x0_p8, p8_x0_names),
        ("p9", p9, x0_p9, p9_x0_names),
        ("p10", p10, x0_p10, p10_x0_names),
        # ("p11", p11, x0_p11,p11_x0_names),
        # ("p12", p12, x0_p12,p12_x0_names),
        ("p11", p13, x0_p13, p13_x0_names),
        ("p12", p14, x0_p14, p14_x0_names),
        ("p13", p15, x0_p15, p15_x0_names),
        ("p14", p16, x0_p16, p16_x0_names),
        ("p15", p17, x0_p17, p17_x0_names),
        ("p16", p18, x0_p18, p18_x0_names),
        # ("p19", p19, x0_p19,p19_x0_names),
        ("p17", p20, x0_p20, p20_x0_names),
        ("p18", p21, x0_p21, p21_x0_names),
        ("p19", p22, x0_p22, p22_x0_names),
        ("p20", p23, x0_p23, p23_x0_names),
        ("p21", p24, x0_p24, p24_x0_names),
        ("p22", p25, x0_p25, p25_x0_names),
        ("p23", p26, x0_p26, p26_x0_names),
        ("p24", p27, x0_p27, p27_x0_names),
        ("p25", p28, x0_p28, p28_x0_names),
        ("p26", p29, x0_p29, p29_x0_names),
        ("p27", p30, x0_p30, p30_x0_names),
        # ("p28", p31, x0_p31,p31_x0_names),
        ("p28", p32, x0_p32, p32_x0_names),
        # ("p33", p33, x0_p33,p33_x0_names),
    ]

end