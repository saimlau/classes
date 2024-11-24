using PlotlyJS

plot(heatmap(
    z=[[1, missing, 30, 50, 1], [20, 1, 60, 80, 30], [30, 60, 1, -10, 20]],
    x=["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"],
    y=["Morning", "Afternoon", "Evening"],
    hoverongaps=false
))