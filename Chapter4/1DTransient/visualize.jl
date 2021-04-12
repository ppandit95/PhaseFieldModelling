using DelimitedFiles,Plots;pyplot()
data = "Initial_Profile.dat"
dfp = @__DIR__
data = string(dfp,"/",data)
D = readdlm(data)
plot(D)
