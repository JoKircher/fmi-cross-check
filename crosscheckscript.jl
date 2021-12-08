using FMI
using CSV, DataFrames
using StatsPlots, Plots

include("cross-check-functions.jl")

# initialize structure
modeltype = ["cs"]
vendors = ["CATIA/R2016x", "Dymola/2019FD01", "DS_FMU_Export_from_Simulink/2.3.0", "FMIToolbox_MATLAB/2.3"]
system = ""

# first entry always time
# MixtureGases, ControlledTemperature, CoupledClutches
# [der_p, der_T], [heatCapacitor_T, switch_Controll], [J1_w,J2_w,J3_w,J4_w]
errorsCATIA = [[0.000005, 0.02, 0.0003], [0.005, 0.06, 1], [0.04, 0.4, 0.2, 0.2, 0.3]]
# MixtureGases, ControlledTemperature, CoupledClutches
# [der_p, der_T], [heatCapacitor_T, switch_Controll], [J1_w,J2_w,J3_w,J4_w]
errorsDymola = [[0.000005, 0.015, 0.0003], [0.005, 0.02, 1], [0.09, 0.4, 0.03, 0.03, 0.08]]
# BouncingBalls_sf, TriggeredSubsystems_sf, TestModel1_sf, TestModel2_sf
# [Ball1_pos,Ball2_pos,Ball1_vel,Ball2_vel], [yRising,yFalling,yEither], [Noise,Noise_Filt,Noise_Filt_Sat,Sine,Sine_Filt,Sine_Filt_Quantized]
errorsDS = [[0.000005, 0.0005, 0.0005, 0.0005, 0.0005], [0.000005, 5e-6, 5e-6, 5e-6], [0.000005, 9, 3.1, 1, 0.00005, 0.00005, 0.11], [0.02, 1, 3]]
# Continuous, IntegrateSignal, EmbeddedCode
# [Out1,Out2[1],Out2[2],Out3,Out4,Out5], [Out1], [Out1]
errorsFMIToolbox = [[0.000005, 1e-20, 1e-20, 1e-20, 1e-15, 1e-16, 1e-14], [0.02, 0.008], [0.000005, 1e-15]]

if Sys.iswindows()
    system = "win$(Sys.WORD_SIZE)"
elseif Sys.islinux()
    system = "lin$(Sys.WORD_SIZE)"
elseif Sys.isapple()
    system = "darwin$(Sys.WORD_SIZE)"
else
    @assert false "OS not supported"
end

for type in modeltype
    for vendor in vendors
        maxerror = []
        if vendor == "CATIA/R2016x" || vendor == "Dymola/2019FD01"
            FMUs = ["MixtureGases", "ControlledTemperature", "CoupledClutches"]
            if vendor == "CATIA/R2016x"
                maxerror = errorsCATIA
            else
                maxerror = errorsDymola
            end
        elseif vendor == "DS_FMU_Export_from_Simulink/2.3.0"
            FMUs = ["BouncingBalls_sf", "TriggeredSubsystems_sf", "TestModel1_sf", "TestModel2_sf"]
            maxerror = errorsDS
        else
            FMUs = ["Continuous", "IntegrateSignal", "EmbeddedCode"]
            maxerror = errorsFMIToolbox
        end
        n = 0
        for fmuname in FMUs
            n = n + 1
            if isfile("fmus/2.0/$type/$system/$vendor/$fmuname/$fmuname.fmu")
                pathtoFMU = "fmus/2.0/$type/$system/$vendor/$fmuname/$fmuname.fmu"
            else
                println("WARNING! FMU not found")
            end

            ref =string(fmuname, "_ref")
            out = string(fmuname, "_out")
            in = string(fmuname, "_in")

            name, rowsres, colsres = readCSV(type, system, vendor, fmuname, ref)

            # initialize simulation
            t_start = rowsres[1].time
            t_stop = rowsres[end].time
            failed = false

            fmu = fmiLoad(pathtoFMU)
            j = 0    
            fmiInstantiate!(fmu;loggingOn=false)
            if isfile("fmus/2.0/$type/$system/$vendor/$fmuname/$in.csv")               
                input, inputrows, inputcols = readCSV(type, system, vendor, fmuname, in)
                if type == "cs"
                    df = fmiSimulatewithInput(fmu, rowsres, input, inputcols, name, t_start, t_stop)
                else
                    @warn "Noch nicht implementiert!"
                end
            else
                if type == "cs"
                    success, data = fmiSimulate(fmu,t_start, t_stop;recordValues=name[2:end], saveat=rowsres.time, setup=true)
                    df = DataFrame()
                    df.time = data.t
                    for i in 2:length(name)
                        insertcols!(df, i, Symbol(name[i]) => [x[i-1] for x in data.saveval])
                    end
                else
                    @warn "Noch nicht implementiert"
                end
            end

            colsdf = eachcol(df)
            rowsdf = eachrow(df)
            
            if isfile("fmus/2.0/$type/$system/$vendor/$fmuname/$in.csv")
                failed, error = calculateErrorwithInput(name, colsdf, colsres, maxerror, fmuname, n)
            else
                failed, error = calculateError(name, colsdf, colsres, maxerror, fmuname, n)
            end
            println(failed)
            if !ispath("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/")
                mkpath("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/")
            end
            if failed
                println("Test failed for $fmuname")
                if isfile("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/passed")
                    rm("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/passed")
                end
                touch("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/failed")
            else
                if isfile("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/failed")
                    rm("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/failed")
                end
                touch("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/passed")
            end
            
            CSV.write("results/2.0/$type/$system/FMI_jl/v0.1.6/$vendor/$fmuname/$out.csv", df)
            fmiUnload(fmu)
        end
    end
end
