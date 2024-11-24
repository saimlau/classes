# Problem 6
using JSON3
using Plots
using LaTeXStrings

# part a
data = JSON3.read("pop_dyn_data.json");

A = copy(data[:A][:data]);
A = stack(A);
n = data[:n][:data];

x = zeros(n);
x[4] = 1;
s = zeros(n);
t = 0;
s[4] = -1;

while any(s.==0.0)
    t +=1;
    x = A*x;
    for i=1:n
        if x[i]!=0.0 && s[i]==0.0
            s[i] = t;
        end
    end
end
s[4] = 0.0;
println("s = $s.T")

# part b
B = A^10;
b1_til = B[1,:];
println("[A^10]_(1st row) = $b1_til")

x_0 = [-1,1,-1,1,1,-1,-1,-1,-1,-1];
temp = B*x_0;
x_1_10 = temp[1];
println("x_1(10) = $x_1_10")

t_s = 0:40;
x_1_s = zeros(length(t_s));
x_1_s[1] = x_0[1]
x = x_0;
for i in t_s[2:end]
    x = A*x;
    x_1_s[i+1] = x[1];
end
plot(t_s,x_1_s);
xlabel!(L"t");
ylabel!(L"x_1(t)");
savefig("p6b.png")
