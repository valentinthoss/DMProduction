dir = @__DIR__

#Import the f(M,a) table from BlackHawk with adjusted quark masses
path_fM_table = string(dir,"/../Data/BlackHawk/fM_nograv.txt")
fM_table = DelimitedFiles.readdlm(path_fM_table)

#Import the g(M,a) table from BlackHawk with adjusted quark masses
path_gM_table = string(dir,"/../Data/BlackHawk/gM_nograv.txt")
gM_table = DelimitedFiles.readdlm(path_gM_table)

#Import the γ(x) table from BlackHawk for various particle spins
path_γ = string(dir,"/../Data/BlackHawk/spin_")
γ_table_0 = DelimitedFiles.readdlm(string(path_γ,"0.txt"))
γ_table_05 = DelimitedFiles.readdlm(string(path_γ,"0.5.txt"))
γ_table_1 = DelimitedFiles.readdlm(string(path_γ,"1.txt"))
γ_table_15 = DelimitedFiles.readdlm(string(path_γ,"1.5.txt"))
γ_table_2 = DelimitedFiles.readdlm(string(path_γ,"2.txt"))
const x_tab = Array{Float64}(γ_table_0[1,2:end])
const γ_0_tab = Array{Float64}(γ_table_0[2:end,2:end])
const γ_05_tab = Array{Float64}(γ_table_05[2:end,2:end])
const γ_1_tab = Array{Float64}(γ_table_1[2:end,2:end])
const γ_15_tab = Array{Float64}(γ_table_15[2:end,2:end])
const γ_2_tab = Array{Float64}(γ_table_2[2:end,2:end])
γ_fits_0 = DelimitedFiles.readdlm(string(path_γ,"0_fits.txt"))
γ_fits_05 = DelimitedFiles.readdlm(string(path_γ,"0.5_fits.txt"))
γ_fits_1 = DelimitedFiles.readdlm(string(path_γ,"1_fits.txt"))
γ_fits_15 = DelimitedFiles.readdlm(string(path_γ,"1.5_fits.txt"))
γ_fits_2 = DelimitedFiles.readdlm(string(path_γ,"2_fits.txt"))
const γ_0_fit = Array{Float64}(γ_fits_0[2:end,2:end])
const γ_05_fit = Array{Float64}(γ_fits_05[2:end,2:end])
const γ_1_fit = Array{Float64}(γ_fits_1[2:end,2:end])
const γ_15_fit = Array{Float64}(γ_fits_15[2:end,2:end])
const γ_2_fit = Array{Float64}(γ_fits_2[2:end,2:end])

const m_t = fM_table[2:end,1]*GeV_to_g #Load masses
const a_t = fM_table[1,2:end] #Load BH spin values
const fM_t = fM_table[2:end,2:end]*GeV_to_g^3*GeV_to_1_s #Load f(M,a)
const gM_t = gM_table[2:end,2:end]*GeV_to_g^3*GeV_to_1_s #Load g(M,a)
const fM_interp = linear_interpolation((m_t,a_t),fM_t,extrapolation_bc=Flat()) #Interpolate
const gM_interp = linear_interpolation((m_t,a_t),gM_t,extrapolation_bc=Flat()) #Interpolate
const m_tab = 10 .^ LinRange(9,18,2000)
a_tab_temp = a_t
for i=1:4
    global a_tab_temp = sort(vcat(a_tab_temp,sqrt.(a_tab_temp[2:end-1].*a_tab_temp[3:end])))
end
const a_tab = a_tab_temp
const fM_tab = fM_interp(m_tab,a_tab);
const gM_tab = gM_interp(m_tab,a_tab);


#HDMSpectra files/data
const i_part = [1,2,3,4,5,6,11,13,15,12,14,16,21,22,23,24,25] #Particle IDs
const s_part = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.0,0.0] #Particle spins
const g_part = [12,12,12,12,12,12,4,2,4,2,4,2,16,2,3,6,1] #Particle d.o.f.
const n_part = length(i_part)
const part = ["Up","Down","Top","Bottom","Strange","Charm","Electron","Muon","Tau","e-Neutrino","m-Neutrino","t-Neutrino","Gluon","Photon","W-Boson","Z-Boson","Higgs"]
const path_tables = string(dir,"/../Packages/HDMSpectra/tables/table.jld2")
const path_hdm = string(dir,"/../Packages/HDMSpectra/data/HDMSpectra.hdf5")

#Degrees of freedom in the early universe
dof_table = readdlm(string(dir,"/../Data/dof_table.csv"),',')
dof_energy_interp = linear_interpolation(reverse(dof_table[:,1]),reverse(dof_table[:,5]),extrapolation_bc=Flat())
dof_entropy_interp = linear_interpolation(reverse(dof_table[:,1]),reverse(dof_table[:,7]),extrapolation_bc=Flat())
a = (dof_table[:,1]*GeV_to_K/2.7255).^(-1) .* (dof_table[:,7]/dof_table[end,7]).^(-1/3)
const rad_correct_x = a
const rad_correct_y = dof_table[:,5]/dof_table[end,5] .* (dof_table[end,7] ./ dof_table[:,7]).^(4/3)