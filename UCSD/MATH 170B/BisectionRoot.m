function Xs = BisectionRoot(Fun,a,b)
if Fun(a)==0, Xs = a; return, 
elseif Fun(b)==0, Xs = b; return,
elseif Fun(a)*Fun(b)>0, error("a and b needs to be on opposite sides of the solution.")
end
if a>b, tem = a; a=b; b=tem; end
epi = 0.000001;
aa = a;
bb = b;
Xs = (aa+bb)/2; 
while abs(Fun(Xs))>=epi
    Xs = (aa+bb)/2;
    switch sign(Fun(Xs)*Fun(aa))
        case 0
            return
        case 1
            aa = Xs;
        case -1
            bb = Xs;
        otherwise
            error("something weird went wrong!!!")
    end
end

end