using DelimitedFiles,Plots;pyplot()
data = ["Initial_Profile.dat" "Output-0.dat" "Output-100.dat" "Output-200.dat" "Output-300.dat" "Output-400.dat" "Output-500.dat" "Output-600.dat" "Output-700.dat" "Output-800.dat"]
dfp = @__DIR__
for i in 1:10
    data[i] = string(dfp,"/",data[i])
    D = readdlm(data[i])
    plot(D)
    savefig(string(chop(data[i],tail=3),"png"))
end
