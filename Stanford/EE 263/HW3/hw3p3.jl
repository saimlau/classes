# Problem 3
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW3/rev_eng_smooth_data.json")
n = data["T"];
u = data["u"];
y = data["y"];

g(x) = [y[x] 2*y[x]-y[x-1]-y[x+1] y[x-2]-4*y[x-1]+6*y[x]-4*y[x+1]+y[x+2]];

i = Int(ceil(n/2))

sol = [g(i-1);g(i);g(i+1)]\[u[i-1]-y[i-1];u[i]-y[i];u[i+1]-y[i+1]]

# guess = [0.1;2;10]

# test1 = [g(5);g(6);g(7)]*sol-[u[5]-y[5];u[6]-y[6];u[7]-y[7]]
# test2 = [g(5);g(6);g(7)]*guess-[u[5]-y[5];u[6]-y[6];u[7]-y[7]]
