# Problem 5
t = collect(1:1000);
z = 5*sin.(t/10 .+ 2) + 0.1 * sin.(t) + 0.1*sin.(2*t .- 5);

z_til = zeros(1000);
z_til[1:3] = z[1:3];

for i=4:1000
    z_til[i] = 3*z_til[i-1]-3*z_til[i-2]+1*z_til[i-3];
end

rms_err = sqrt(sum((z_til[4:1000].-z[4:1000]).^2)./sum(z[4:1000].^2))
