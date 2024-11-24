# Problem 1
using LinearAlgebra
using Plots
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW5/tomo_data.json")
N = data["N"]
n = data["npixels"]
y = data["y"]
L = data["line_pixel_lengths"]

x = transpose(L)\y;
A = reshape(x,n,n)
heatmap(A, yflip=true, aspect_ratio=:equal, color=:gist_gray,
cbar=:none, framestyle=:none)
savefig("HW5/p1.png")