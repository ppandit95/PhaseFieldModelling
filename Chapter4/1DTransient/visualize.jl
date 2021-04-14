using DelimitedFiles,Plots;pyplot()
data = ["Initial_Profile.dat" "Output-0.dat" "Output-1.dat" "Output-2.dat" "Output-3.dat" "Output-4.dat" "Output-5.dat" "Output-6.dat"]
dfp = @__DIR__
data = string(dfp,"/",data)
D = readdlm(data)
plot(D)
savefig(string(chop(data,tail=3),"png"))
