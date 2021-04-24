using DelimitedFiles,Plots;pyplot()
data0 = readdlm("/media/pushkar/Data/Personal Files/PhaseFieldModelling/Chapter4/CaseStudyI/Initial_Profile.dat")
#=data1 = readdlm("/home/pushkar/Task_assignment_pics/AllenCahn2DFDM/Output-10.dat")
data2 = readdlm("/home/pushkar/Task_assignment_pics/AllenCahn2DFDM/Output-20.dat")
data5 = readdlm("/home/pushkar/Task_assignment_pics/AllenCahn2DFDM/Output-50.dat")
data10 = readdlm("/home/pushkar/Task_assignment_pics/AllenCahn2DFDM/Output-100.dat")=#
heatmap(data0,cgrad=[:blue :red])
savefig("/media/pushkar/Data/Personal Files/PhaseFieldModelling/Chapter4/CaseStudyI/Initial_Profile.png")
#=heatmap(data1,c=:greys)
savefig("/home/pushkar/Task_assignment_pics/AllenCahn2DFDM/Output-10.png")
heatmap(data2,c=:greys)
savefig("/home/pushkar/Task_assignment_pics/AllenCahn2DFDM/Output-20.png")=#
