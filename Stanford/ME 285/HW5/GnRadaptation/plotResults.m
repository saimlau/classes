clear, clc, close all
names = ["4P", "1-2Q", "4P_1-2Q"];
for name=names
    data = readmatrix("simulation_"+name+"_output");
    t = data(:,1);
    essure = data(:,2);
    flow = data(:,3);
    a = data(:,4);
    h = data(:,5);
    
    figure(1)
    subplot(1,2,1)
    plot(t,a)
    title("Radius")
    xlabel("t (days)")
    subplot(1,2,2)
    plot(t,h)
    title("Thickness")
    xlabel("t (days)")
    saveas(gcf, name+".png")
end 


