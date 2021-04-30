using DelimitedFiles,Plots;pyplot()
DIR = @__DIR__
for i in 0:1000:9000
    path = string(DIR,"/Concentration_Profile_t_",i,".dat")
    data = reshape(readdlm(path)[:,3],64,64)
    heatmap(data,cgrad=[:blue :red])
    savefig(string(DIR,"/Results/Concentration_Profile_t_",i,".png"))
end
 path = string(DIR,"/Gibbs_Energy.dat")
 data = readdlm(path)
 plot(1:length(data),data,xlabel="Time",ylabel="Total Gibbs Energy",legend=false)
 savefig(string(DIR,"/Results/Gibbs_Energy.png"))
