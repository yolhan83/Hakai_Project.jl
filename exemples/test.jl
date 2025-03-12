using Hakai_Project

@profview try
    Hakai_Project.main("./exemples/input/car-crash-N2k.inp")
catch e
    println(e)
end
