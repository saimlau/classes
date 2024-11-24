using LinearAlgebra

G = [100,10,-1]
G_5 = G./norm(G)*5

th = [-1,1,1]
th_new = th+0.5*G_5