using FMI
using CSV, DataFrames

function readCSV(type, system, vendor, fmuname, ref)
    if isfile("fmus/2.0/$type/$system/$vendor/$fmuname/$ref.csv")
        df = CSV.read("fmus/2.0/$type/$system/$vendor/$fmuname/$ref.csv", DataFrame)
    else
        println("WARNING! Reference-CSV not found!")
    end

    # get column name
    name = names(df)

    # get iterators for rows and columns
    rowsdf = eachrow(df)
    colsdf = eachcol(df)
    [name, rowsdf, colsdf, df]
end

function fmiSimulatewithInput(fmu, refrows, input, inputcols, name, t_start, t_stop)
    # get the initial input values
    inputindex = 1
    inputatpoint = Array{FMI.fmi2Real}(undef, 0)
    for i in 1:length(inputcols)
        push!(inputatpoint, inputcols[i][inputindex])
    end

    # initialize FMU
    fmi2EnterInitializationMode(fmu)
    fmiSetReal(fmu,input[2:end], inputatpoint[2:end])
    fmi2ExitInitializationMode(fmu)

    fmiSetupExperiment(fmu, t_start)

    # select csv with smaller timesteps for simulation time steps
    isinput = false
    if fmiGetModelName(fmu) == "IntegerNetwork1"
        if refrows[4].time - refrows[3].time > inputcols[1][4]-inputcols[1][3]
            isinput = true
        end
    else
        if refrows[2].time - refrows[1].time > inputcols[1][2]-inputcols[1][1]
            isinput = true
        end
    end

    saveat = []
    if fmiGetModelName(fmu) == "IntegerNetwork1"
        dt = min((refrows[4].time - refrows[3].time),(inputcols[1][4]-inputcols[1][3]))
    else
        dt = min((refrows[2].time - refrows[1].time),(inputcols[1][2]-inputcols[1][1]))
    end
    t = t_start
    if isinput
        for i in 1:length(refrows)
            push!(saveat, refrows[i].time)
        end
    else
        step = inputcols[1][2] - inputcols[1][1]
        for i in 1:length(inputcols[1])
            push!(saveat, i * step)
        end
    end
    # create dataframe for storage
    help = zeros(2, length(name))
    sim = DataFrame(help, :auto)
    rename!(sim, map(Symbol, name))
    delete!(sim, 1)
    delete!(sim, 1)

    # value represents the output values at the current time
    value = Array{FMI.fmi2Real}(undef, length(name)-1)
    saveatidx = 1
    while t < t_stop
        fmiDoStep(fmu, t, dt)
        if isinput
            # add current values to simulation result data frame (only works if in and ref csv time steps are multiples of each other)
            if abs(t - saveat[saveatidx]) < 1e-4
                if length(name) == 2
                    value = [t, fmiGetReal(fmu, name[2])]
                else
                    if fmiGetModelName(fmu) == "TestModel2_sf"
                        push!(value, t)
                        push!(value, fmiGetBoolean(fmu, name[2]))
                        push!(value, fmiGetReal(fmu, name[3]))
                    elseif fmiGetModelName(fmu) == "IntegerNetwork1"
                        value = fmiGetInteger(fmu, name[2:end])
                        pushfirst!(value, t)
                    else
                        value = fmiGetReal(fmu, name[2:end])
                        pushfirst!(value, t)
                    end
                end
                push!(sim, Tuple(value))
                if saveatidx != length(saveat)
                    saveatidx += 1
                end
            end
            if fmiGetModelName(fmu) == "IntegerNetwork1"
                fmiSetInteger(fmu,input[2:end], Array{FMI.fmi2Integer}(inputatpoint[2:end]))
            else
                fmiSetReal(fmu, input[2:end], inputatpoint[2:end])
            end
            if inputindex < length(inputcols[1])
                inputindex = inputindex + 1
                # update inputs for next iteration
                for i in 1:length(inputcols)            
                    inputatpoint[i] = inputcols[i][inputindex]
                end
            end
        else
            # for variables in fmu.modelDescription.modelVariables
            #   println(variables.datatype.datatype)
            # end
            if length(name) == 2
                value = [t, fmiGetReal(fmu, name[2])]
            else
               if fmiGetModelName(fmu) == "TestModel2_sf"
                    value = [t]
                    push!(value, fmiGetBoolean(fmu, name[2]))
                    push!(value, fmiGetReal(fmu, name[3]))
                elseif fmiGetModelName(fmu) == "IntegerNetwork1"
                    value = fmiGetInteger(fmu, name[2:end])
                    value = convert(Array{Float64}, value)
                    pushfirst!(value, t)
                else
                    value = fmiGetReal(fmu, name[2:end])
                    pushfirst!(value, t)
                end
            end
            push!(sim, Tuple(value))
            if saveatidx <= length(saveat)
                while abs(saveat[saveatidx] - inputcols[1][inputindex]) > 1e-5 && inputindex < length(inputcols[1])
                    inputindex += 1
                end
                
            end
            if abs(t - saveat[saveatidx]) < 1e-5
                if fmiGetModelName(fmu) == "IntegerNetwork1"
                    fmiSetInteger(fmu,input[2:end], Array{FMI.fmi2Integer}(inputatpoint[2:end]))
                else
                    fmiSetReal(fmu, input[2:end], inputatpoint[2:end])
                end
                if inputindex < length(inputcols[1])
                    inputindex = inputindex + 1
                     # update inputs for next iteration
                    for i in 1:length(inputcols)            
                        inputatpoint[i] = inputcols[i][inputindex]
                    end
                end
                if saveatidx != length(saveat)
                    saveatidx += 1
                end
            end
        end
        t += dt
    end
    sim
end

function calculateErrorwithInput(name, colssim, colsref, maxerror, fmuname, n)
    error = Matrix{Float64}(undef, length(name), length(colssim[1]))
    failed = false
    #calculate errors for all model variables
    for i in 1:length(colssim)
        k = 2
        for j in 1:length(colssim[1])
            if k <= length(colsref[1]) 
                if (colsref[1][k-1] == colsref[1][k])
                    idx = k
                    while colsref[1][k] == colsref[1][idx] && idx < length(colsref[1])
                        idx = idx + 1
                    end
                    error[i, j] = min(abs(colssim[i][j] - colsref[i][k-1]), abs(colssim[i][j]- colsref[i][idx]))
                    k = idx + 1
                else
                    error[i,j] = abs(colssim[i][j] - colsref[i][k-1])
                    k = k + 1
                end
            else
                error[i,j] = abs(colssim[i][j] - colsref[i][k-1])
            end
            
            if error[i,j] > maxerror[n][i]
                println(name[i])
                failed = true
            end
            
        end
    end
    failed, error
end

function calculateError(name, colssim, colsref, maxerror, fmuname, n)
    error = Matrix{Float64}(undef, length(name), length(colssim[1]))
    failed = false
    #calculate errors for all model variables
    for i in 1:length(colssim)
        k = 2
        for j in 1:length(colssim[1])
            if k <= length(colssim[1]) 
                if (colssim[1][k-1] == colssim[1][k])
                    idx = k
                    while colssim[1][k-1] == colssim[1][idx] && idx < length(colssim[1])
                        idx = idx + 1
                    end
                    error[i, j] = min(abs(colssim[i][k-1] - colsref[i][k-1]), abs(colssim[i][k-1]- colsref[i][idx]), abs(colssim[i][idx] - colsref[i][k-1]), abs(colssim[i][idx]- colsref[i][idx]))
                    k = idx + 1
                else
                    error[i,j] = abs(colssim[i][k] - colsref[i][k])
                    k = k + 1
                end
            else
                error[i,j] = abs(colssim[i][k-1] - colsref[i][k-1])
            end
            
            if error[i,j] > maxerror[n][i]
                failed = true
            end
            
        end
    end
    failed, error
end
