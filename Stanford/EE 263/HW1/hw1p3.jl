# Promblem 3
using JSON3

data = JSON3.read("ts_data.json");
A = copy(data[:A][:data]);
A = stack(A);
K = data[:K][:data];  # Number of time slots
n = data[:n][:data];  # Number of nodes

# Part a
function create_e(active_node)
    x = zeros(n);
    x[active_node] = 1;
    return x
end
function update(x, count)
    tem = count%3;
    tem==0 ? tem=1 : tem=tem;
    B = A.==tem;
    return B*x
end
function find_path(x, path, count, n_goal)
    x_new = update(x, count)
    if x_new[n_goal]==1
        append!(path, n_goal)
        println("message arrived at $n_goal,")
        println("Path: $path")
        return path
    else
        indx = findall(y->y==1, x_new)
        for i in indx
            tem = create_e(i);
            new_path = copy(path);
            push!(new_path, i)
            return find_path(tem, new_path, count+1, n_goal)
        end
    end
end
function find_path(n_init, n_goal)
    x = create_e(n_init)
    path_init = [n_init];
    return find_path(x, path_init, 1, n_goal)
end

instruction = find_paths(5, 18)

# # while x[n_goal]==0
# for i =1:1000
#     for k=1:K
#         x_backup = x;
#         t += 1;
#         B = A.==k;
#         x = B*x;
#         if all(x.==0.0)
#             x = x_backup;
#         end
#     end
# end
# println(t)

# Part b
x = zeros(n); # States
x[7] = 1; # init message
t = 0;
x_last = x;
while any(x.==0.0)
    for k=1:K
        t += 1;
        tem = (A.==k)*x
        x += tem;
        x = (x.≠0.0)
        println(x)
    end
    if x_last==x
        t -= K;
        break
    end
    x_last = x;
end
println(x)
println(t)


function testing(stop_factor = 3)
    println("Getting $stop_factor")
    if stop_factor<0
        println("Done")
        return true
    end
    testing(stop_factor-1)
    testing(-6)
    return false
end
testing()
