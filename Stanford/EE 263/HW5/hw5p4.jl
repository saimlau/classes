# Problem 4
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW5/tempdata.json")
x_train = data["x_train"]
y_train = data["y_train"]
x_test = data["x_test"]
y_test = data["y_test"]

n = length(y_train)
scatter(x_train,y_train, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12, label="data points",
background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom)
xlabel!("x")
ylabel!("y")

function f(j,x,m)
    if x<(j-1)*12/m
        return 0
    elseif x<j*12/m
        return m*x/12-(j-1)
    else
        return 1
    end
end

xs = 0:0.01:12
ms = [3,6,12,24]
err = zeros(length(ms))
test_err = zeros(length(ms))
for r in eachindex(ms)
    m = ms[r]
    F = ones((n,m+1))
    for i in 1:n
        for j in 2:m+1
            xi = x_train[i]
            jj = j-1
            F[i,j] = f(jj,xi,m)
        end
    end
    a = F\y_train
    g(x) = sum([a[j]*f(j-1,x,m) for j in 1:m+1])
    plot!(xs, g.(xs), label="m=$m", linewidth = 2)
    err[r] = norm(y_train-F*a)^2
    test_err[r] = norm(y_test-g.(x_test))^2
end
savefig("HW5/p4c.png")

plot(ms,err, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12, label="squared 2-norm error",
background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, linewidth = 2)
xlabel!("m")
ylabel!("error")
savefig("HW5/p4d.png")

plot(ms,test_err, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12, label="test error",
background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, linewidth = 2)
xlabel!(L"m")
ylabel!(L"J^{test}")
savefig("HW5/p4e.png")


