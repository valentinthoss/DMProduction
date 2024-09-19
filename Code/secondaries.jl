#Import secondaries from BlackHawk and HDMSpectra
m_space = 10 .^ LinRange(2,19,150) #Set mass space for contraints
m_path = string(dir,"/../Packages/BlackHawk/version_finale/m_space.txt") #Set path for mass data
emin = 1e-8 #Minimum energy in GeV
emax = 1e15 #Maximum energy in GeV
enum = 2000 #Number of energy values for primary spectrum
DelimitedFiles.writedlm(m_path,m_space) #Save mass space to disk
dir = pwd()
cd(string(dir,"/../Packages/BlackHawk/version_finale"))
cmd = `./run.sh m_space.txt $emin $emax $enum`
if compute_spectra
    run(cmd)
end
cd(dir)

const hadro = ["Pythia_BBN","Herwig_BBN","Pythia","Hazma","HDMSpectra"]

path = string(dir,"/../Packages/BlackHawk/version_finale/results/")
primary_file = string(path,"primary/primary_1.txt")
tmp,header_prim = DelimitedFiles.readdlm(primary_file,skipstart=1,header=true)
header_prim = header_prim
n_prim = size(tmp)[2]

primary_data = zeros(length(m_space),n_prim,enum)
for i=1:length(m_space)
    primary_file = string(path,"primary/primary_",i,".txt")
    primary_data[i,:,:] .= transpose(DelimitedFiles.readdlm(primary_file,skipstart=2))
end
    
secondary_data = []
for j=1:length(hadro)
    secondary_file = string(path,hadro[j],"/secondary_1.txt")
    if hadro[j] !== "HDMSpectra"
        tmp = DelimitedFiles.readdlm(secondary_file,skipstart=2)
    else
        tmp = DelimitedFiles.readdlm(secondary_file)
    end
    push!(secondary_data,zeros(length(m_space),size(tmp)[2],size(tmp)[1]))
    for i=1:length(m_space)
        secondary_file = string(path,hadro[j],"/secondary_",i,".txt")
        if hadro[j] !== "HDMSpectra"
            secondary_data[end][i,:,:] .= transpose(DelimitedFiles.readdlm(secondary_file,skipstart=2))
        else
            secondary_data[end][i,:,:] .= transpose(DelimitedFiles.readdlm(secondary_file))
        end
    end
end

const int_γ = [generate_spectrum_interpolation(secondary_data[i],m_space,col=2) for i=1:length(secondary_data)]
const int_e = [generate_spectrum_interpolation(secondary_data[i],m_space,col=3) for i=1:length(secondary_data)]
    
function spectrum_particle_e_sec(e,m;h="select",part_fin="Photon")
    if part_fin == "Photon"
        interp = int_γ
        s = 1.0
        g = 2.0
    elseif part_fin == "Electron"
        interp = int_e
        s = 0.5
        g = 4.0
    else
        println("Particle not found")
    end
    if h=="select"
        if m < 1.057e11
            j = 5
        elseif m < 1.057e14
            j = 2
        else
            j = 4
        end
        return spectrum_interpolate(interp[j],m,e,s=s,g=g)/g
    else
        j = findfirst(x->x==h,hadro)
        return spectrum_interpolate(interp[j],m,e,s=s,g=g)/g
    end
end