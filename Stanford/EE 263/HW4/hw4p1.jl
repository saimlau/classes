# Problem 1
using LinearAlgebra
using Plots
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW4/leave_one_out_cross_validation_data.json")
N = Int(data["N"])
t = data["t"]
z = data["z"]

param = []
error = zeros(9)
LOOCV_error = zeros(9)

for d in 1:9
    A = ones((N,d+1))
    for i in 1:d
        A[:,i+1] = t.^i
    end
    x = A\z
    H = A*inv(transpose(A)*A)*transpose(A)
    Hii = diag(H)
    push!(param,x)
    error[d] = norm(z-A*x)^2
    LOOCV_error[d] = norm((z-A*x)./(Hii.-1))^2
    if d==2
        L_err_simp = LOOCV_error[d]
        L_err = 0.0
        for i in 1:N
            A_til = A[(1:N).!=i,:]
            y_i_til = z[(1:N).!=i]
            x_i_til = A_til\y_i_til
            r_i_til = z[i]-dot(A[i,:],x_i_til)
            L_err += r_i_til^2
        end
        println("For d=2,")
        println("Direct way, L_err = $L_err")
        println("Inirect way, L_err = $L_err_simp")
    end
end

plot(1:9, error, label="squared fitting err")
plot!(1:9, LOOCV_error, label="LOOCV err")
xlabel!("d")
savefig("HW4/err-d.png")

err4 = error[4]
LOOCV_err4 = LOOCV_error[4]
println("For d=4,")
println("Error = $err4")
println("LOOCV Error = $LOOCV_err4")
scatter(t,z,label="data")
ts = minimum(t)-0.3:0.01:maximum(t)+0.3
pred = [evalpoly(ts[i],param[4]) for i in eachindex(ts)]
plot!(ts,pred, label="d=4 fit")
xlabel!("t")
ylabel!("z")
savefig("HW4/d4_fit.png")

