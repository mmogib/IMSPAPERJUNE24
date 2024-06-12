# # Example usage:
# folder = "paperCode/data"
# files = map(x -> "$folder/$x", readdir(folder))
# # fff = map(x -> begin
# #         sp = splitext("$(folder)/$(basename(x))")
# #         "$(sp[1])_trimmed$(sp[2])"
# #     end, files)
# # foreach(enumerate(files)) do (i, f)
# #     lns = strip.(readlines(f))
# #     open(fff[i], "w") do io
# #         writedlm(io, lns)
# #     end
# # end
# T = Array{Float64,2}(undef, 208, 6)
# for (i, f) in enumerate(files)
#     d1 = readdlm(f, Float64, header=true)
#     bb1 = d1[1][:, 1] .!= 0
#     d1[1][bb1, 2:3] .= NaN
#     T[:, i] = d1[1][:, 3]
#     println(i, f)
# end
