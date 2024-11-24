# Problem 1
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW6/gauss_fit_data.json");
N = data["N"]
t = data["t"]
y = data["y"]

ee(pp) = exp.(-(t.-pp[2]).^2.0./pp[3].^2);  # Saving some typing
r(p) = p[1].*ee(p).-y;

function GaussFit(init_p)
    fig = scatter(t,y, label="data", layout=(1,2), size = (660, 360))
    p = init_p;  # Initial guess, [a,mu,sigma]
    plot!(fig, t,p[1].*ee(p), label="Initial Guess", seriestype=:path, linestyle=:dash, linewidth=2, color = "red", subplot=1)
    title!(fig, "Initial [a,μ,σ] = $p", subplot=1)

    max_iter = 200;
    Dr = zeros(N,3);
    E = [];
    for k in 1:max_iter
        tem = ee(p)
        Dr[:,1] = tem;
        Dr[:,2] = 2.0.*p[1].*(t.-p[2])./p[3].^2.0.*tem;
        Dr[:,3] = 2.0.*p[1].*(t.-p[2]).^2.0./p[3].^3.0.*tem;
        p = inv(transpose(Dr)*Dr)*transpose(Dr)*(Dr*p-r(p))
        push!(E, norm(r(p))/sqrt(N))
        if k>=2 && abs(E[end]-E[end-1])<1e-5
            plot!(fig, 1:k, E, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
            linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom, color = "green", subplot=2)
            xlabel!(fig, "k", subplot=2)
            ylabel!(fig, "E", subplot=2)
            title!(fig, "RMS error", subplot=2)
            break
        end
    end
    plot!(fig, t,p[1].*ee(p), label="Result", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom, color = "green", subplot=1)
    xlabel!(fig, "t", subplot=1)
    ylabel!(fig, "y", subplot=1)
    savefig("HW6/$init_p.png")
end

GaussFit([12,50,25]);
GaussFit([15,60,20]);
GaussFit([3,90,1]);

