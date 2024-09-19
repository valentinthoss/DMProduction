dir = @__DIR__

path_flux = string(dir,"/../Data/Flux data/")

#Loading Isatis data
flux_gc_files = [
"flux_INTEGRAL_1107.0200.txt",
"flux_COMPTEL_1107.0200.txt",
"flux_EGRET_9811211.txt",
"flux_Fermi-LAT_1101.1381.txt"
]

flux_gc_files = string.("Isatis/",flux_gc_files)

flux_gc_labels = [
"INTEGRAL",
"COMPTEL",
"EGRET",
"FERMI-LAT",
"LHAASO"
]

flux_exgb_files = [
"flux_HEAO_balloon_9903492.txt",
"flux_COMPTEL_1502.06116.txt",
"flux_EGRET_0405441.txt",
"flux_Fermi-LAT_1410.3696.txt"
]

flux_exgb_files = string.("Isatis/",flux_exgb_files)

flux_exgb_labels = [
"HEAO+baloon",
"COMPTEL",
"EGRET",
"FERMI-LAT",
"LHAASO"
]

flux_gc_isatis = [
    DelimitedFiles.readdlm(string(path_flux,flux_gc_files[i]),skipstart=4,header=true)[1]
    for i=1:length(flux_gc_files)
]

flux_exgb_isatis = [
    DelimitedFiles.readdlm(string(path_flux,flux_exgb_files[i]),skipstart=4,header=true)[1]
    for i=1:length(flux_gc_files)
]

#Load Lhaaso data

flux_gc_lhaaso = DelimitedFiles.readdlm(string(path_flux,"lhaaso.csv"),',')
const flux_gc_data = vcat(flux_gc_isatis,[flux_gc_lhaaso])

flux_exgb_lhaaso = DelimitedFiles.readdlm(string(path_flux,"lhaaso_exgb.csv"),',')
const flux_exgb_data = vcat(flux_exgb_isatis,[flux_exgb_lhaaso]);

const flux_gc_l = [
[-30,30],
[-30,30],
[-30,30],
[-30,30],
[15,125]
]*deg_to_rad

const flux_gc_b = [
[-15,15],
[-15,15],
[-5,5],
[-10,10],
[-5,5]
]*deg_to_rad