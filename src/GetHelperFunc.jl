###############################
#### Independent functions ####
###############################
//
    """
        even(N::Union{Float64,Int64})

        Take any value N and returns the nearest even number. 
    """
    even(N::Union{Float64,Int64}) = 2*Int(div(N,2,RoundNearest))

    Base.:*(par::Parameters,scalar::Real) = Parameters(par.λ₀ * scalar, par.λ₁ * scalar, par.h * scalar)
    Base.:+(par::Parameters,scalar::Real) = Parameters(par.λ₀ + scalar, par.λ₁ + scalar, par.h + scalar)
    Base.:*(dis::DistanceExperSampleStruct,scalar::Real) = 
        DistanceExperSampleStruct(round.(dis.experimental .* scalar,digits=2),round.(dis.sample .* scalar,digits=2))
    

    
    """
    t⃗(var::Variables)

    Takes the var.T and var.Δt to return a StepRangeLen.
    """
    t⃗(var::Variables) = range(0,var.T,step=var.Δt)
    
    """
    x⃗(var::Variables)

    Takes the var.X and var.Δx to returns a StepRangeLen.
    """
    x⃗(var::Variables) = range(0,var.X,step=var.Δx)


    

    """
    Base.range(var::Variables)

    Generate range for spatial domain with 0 and centre. Distretised by Δx.
    """
    

    Base.range(var::Variables) = range(-var.X/2,var.X/2,step=var.Δx)
    Base.range(X::Union{Int64,Float64},var::Variables) = range(-X/2,X/2,step=var.Δx)
    
    
    """
        dimensions(var::Variables)

        Take X (width) and aspect to return a NamedTuple.
    """
    function dimensions(var::Variables)
        if var.Δx != 1
            return Domain(var.X,var.X*var.aspect)
        else
            return Domain(var.X,even(var.X*var.aspect))
        end
    end
    
 
    """
        Δ(par::Parameters,var::Variables)

        Take the parameters (Parameters) and the variables (Variables) to return 
        there discretisation value.

    """
    function Δ(par::Parameters,var::Variables)
        return par*var.Δt
    end

    

    """
        Δ(io,Δ)

        Take an input (io) and a discretisation (Δ) to return its ioΔ value.
    """
    function Δ(io,Δv)
        return io*Δv
    end



    """
        s₀()    

        Draw a single random 0 or 1.
    """
    s₀() = rand((0,1))


    """
        LitRate(Diffusion,var::Variables)

        Take a diffusion coefficient found in literature and convert to the hopping rate.
    """
    function LitRate(Diffusion,var::Variables)
        return Diffusion/var.Δx^2
    end

    """
    LitRate(Diffusion,Δx)

    Take a diffusion coefficient found in literature and convert to the hopping rate.
    """
    function LitRate(Diffusion,Δx)
        return Diffusion/Δx^2
    end


    function param_depend_distplacement(par::Parameters,var::Variables)
        return √((par.h*var.Δx^2)/(1 + (par.λ₁/par.λ₀)))
    end
    

    """
    Get value of parameter or variable in relations to the property of the stuct.
    """
    function get_value(data,data_type::String, data_value::String)
            return @chain data begin
                            getproperty(_,Symbol(data_type))
                            getproperty(_,Symbol(data_value))
                    end
    end         


    """
    Takes type `Vector` of data which is type `SolutionVarParDom` and extracts the mean of the sample and experimental per each parameters
    """
    function get_mean_distance(data,iter,data_type, data_value)
            return @chain data[iter] begin
                            map(x -> read_solution_from_memory(x,SolutionVarParDom),_)
                            map(x -> (get_value(x,data_type, data_value),mean(x.experimental),mean(x.sample)),_)
                            sort(_)
                    end

    end

    """
    Convert output from `load_mean_data` into a `NamedTuple` of vectors.
    """
    function get_mean_vector(data,iter)        
            param = [i[1] for i in data[iter]]
            exper = [i[2] for i in data[iter]]
            sampl = [i[3] for i in data[iter]]
            return (parameter = param, experimental = exper, sample = sampl)        
    end

    """
    Convert output from `load_mean_data` into a `NamedTuple` of vectors, as well the theoretical
    """
    function get_solutions_vector(data,theoretical,iter)        
            param = [i[1] for i in data[iter]]
            exper = [i[2] for i in data[iter]]
            sampl = [i[3] for i in data[iter]]
            theor = [i for i in theoretical[iter]]
            return (parameter = param, experimental = exper, sample = sampl,theoretical = theor)        
    end

    """
    Convert output from `SolutionVarParDom` into a `NamedTuple` of vectors, as well the theoretical
    """
    function get_solutions_vector(data,theoretical)        
            param = [i[1] for i in data]
            exper = [i[2] for i in data]
            sampl = [i[3] for i in data]
            theor = [i for i in theoretical]
            return (parameter = param, experimental = exper, sample = sampl,theoretical = theor)        
    end



    """
    Apply * to the Domain struct
    """
    Base.:*(d::Domain,x) = Domain(d.width * x, d.length * x)

    """ 
        Apply * to T6ssDistance
    """
    Base.:*(d::T6ssDistance,x) = T6ssDistance(d.distance * x)


    """
        Apply * to Walker2D
    """
    Base.:*(w::Walker2D,x) = Walker2D(w.x*x,w.y*x)

    """
        Apply * to WalkerFirstLast
    """
    Base.:*(wfl::WalkerFirstLast,x) = WalkerFirstLast(wfl.w₁ * x,wfl.w₂ * x)
    
//    


