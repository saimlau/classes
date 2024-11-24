using CSV
using DataFrames
using Plots
using Statistics
print("\033c")
dataf = CSV.read("testingRCR/testingRCR-converted-results_1d/testingRCR_flow.csv", DataFrame)
datap = CSV.read("testingRCR/testingRCR-converted-results_1d/testingRCR_pressure.csv", DataFrame)

tf = dataf[:,1];
tp = datap[:,1];
n = length(dataf[1,:])
m = length(tf)
inlet = "branch0_seg0"
outlets = ["branch9_seg0", "branch28_seg0", "branch22_seg0", "branch2_seg0", "branch4_seg2", "branch15_seg0", "branch25_seg0", "branch6_seg0", "branch7_seg0", "branch3_seg0", "branch26_seg0", "branch17_seg2", "branch20_seg0", "branch5_seg0", "branch11_seg2", "branch21_seg0", "branch16_seg0", "branch23_seg0", "branch27_seg0"]
# for seg_to_plot in push!(copy(outlets),inlet)
#     # seg_to_plot = 9
#     plot(tf[m-84:end], dataf[!,seg_to_plot][m-84:end].*60/1000, framestyle = :box, xguidefontsize=15, yguidefontsize=15,legendfontsize=15, ytickfontsize = 15, xtickfontsize = 15, linewidth = 2, label=split(seg_to_plot,"_")[1])
#     xlabel!("time sec")
#     ylabel!("Flow (L/min)")
#     savefig("1DSim1-converted-results_1d/$(seg_to_plot)"*"_flow")

#     plot(tp[m-84:end], datap[!,seg_to_plot][m-84:end]./1333.2, framestyle = :box, xguidefontsize=15, yguidefontsize=15,legendfontsize=15, ytickfontsize = 15, xtickfontsize = 15, linewidth = 2, label=split(seg_to_plot,"_")[1])
#     xlabel!("time (sec)")
#     ylabel!("Pressure (mmHg)")
#     savefig("1DSim1-converted-results_1d/$(seg_to_plot)"*"_pressure")
#     println("Averaged $seg_to_plot pressure ≈ $(mean(datap[!,seg_to_plot])/1333.2) mmHg, ($(mean(datap[!,seg_to_plot])) dyne/cm^2)")
# end

# println("Averaged inlet pressure ≈ $(mean(datap[!,inlet])/1333.2) mmHg, ($(mean(datap[!,inlet])) dyne/cm^2)")
# println("Systolic inlet pressure ≈ $(maximum(datap[!,inlet])/1333.2) mmHg, ($(maximum(datap[!,inlet])) dyne/cm^2)")
# println("Diastolic inlet pressure ≈ $(mean(datap[!,inlet])/1333.2) mmHg, ($(mean(datap[!,inlet])) dyne/cm^2)")

plot()
for seg_to_plot in names(dataf)[2:end]
# for seg_to_plot in outlets
    plot!(tf[m-84:end], dataf[!,seg_to_plot][m-84:end].*60/1000, framestyle = :box, xguidefontsize=15, yguidefontsize=15,legendfontsize=15, ytickfontsize = 15, xtickfontsize = 15, linewidth = 2, label=split(seg_to_plot,"_")[1])
    xlabel!("time sec")
    ylabel!("Flow (L/min)")
end
savefig("testingRCR/testingRCR-converted-results_1d/result"*"_flow")

plot()
for seg_to_plot in names(datap)[2:end]
# for seg_to_plot in outlets
    plot!(tp[m-84:end], datap[!,seg_to_plot][m-84:end]./1333.2, framestyle = :box, xguidefontsize=15, yguidefontsize=15,legendfontsize=15, ytickfontsize = 15, xtickfontsize = 15, linewidth = 2, label=split(seg_to_plot,"_")[1])
    xlabel!("time (sec)")
    ylabel!("Pressure (mmHg)")
end
savefig("testingRCR/testingRCR-converted-results_1d/result"*"_pressure")
