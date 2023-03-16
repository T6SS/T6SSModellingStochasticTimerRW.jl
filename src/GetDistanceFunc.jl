#########################################
#### Distance and Boundary functions ####e 
#########################################
//
 """
    Take domain, amount of iterations and variables for a given simulation and generate a stochastic sample distance distribution.
    Distance distribution will be conditional on the boundary conditions.
 """
 
    function get_sample_distribution(var::Variables,dom::Domain,iter::Iterations)
            RandW() = rand(Int(-(1/2)*dom.width):Int((1/2)*dom.width),10^iter.simulation_iterations)
            RandL() = rand(Int(-(1/2)*dom.length):Int((1/2)*dom.length),10^iter.simulation_iterations)

            if var.boundary == "no_boundary"
                X = colwise(Cityblock(), RandW(),RandW())
                Y = colwise(Cityblock(), RandL(),RandL())
            elseif var.boundary == "double_periodic"
                X = colwise(Cityblock(), RandW(),RandW())
                X = [(x > dom.width / 2 ? dom.width - x : x) for x in X]
                Y = colwise(Cityblock(), RandL(),RandL())
                Y = [(y > dom.length / 2 ? dom.length - y : y) for y in Y]
            elseif var.boundary == "no_flux"
                X = colwise(Cityblock(), RandW(),RandW())
                X = [(x > dom.width / 2 ? dom.width - x : x) for x in X]
                Y = colwise(Cityblock(), RandL(),RandL())
            elseif var.boundary == "cell_boundary"
                X = colwise(Cityblock(), RandW(),RandW())
                X = [(x > dom.width / 2 ? dom.width - x : x) for x in X]
                Y = colwise(Cityblock(), RandL(),RandL())
            else
                #nothing
            end
            return X .+ Y
    end
    
    """
    Get the distance if the boundary is periodic. This shifts the distance by 1/2 the domain in the direction of such distance. This function will accept the distance and domain in either an Int64 or a Float64 as a union.
    """ 
    function get_periodic_boundary_distance(distance::T,domain::T) where T <: Union{Int64,Float64}
        return distance > domain/2 ? domain - distance : distance
    end
    """
    Takes the distance and the boundary, but returns the distance as is.  This function will accept the distance and domain in either an Int64 or a Float64 as a union.
    """
    function get_no_periodic_boundary_distance(distance::T,domain::T) where T <: Union{Int64,Float64}
        return distance
    end
    

    """
    If only one walker with type Walker2D, then return the summed absolute value of the walker. Takes a Walker2D type as input and returns the same type as the walker input.
    """
    function get_distance(w::Walker2D)
        return abs(w.x) + abs(w.y)
    end

    """
    Takes WalkerFirstLast type, which is a pair of Walker2D, say (W,W) -> D, where D is of type Distnace, which is a tuple of the x and y distances each.
    """
    function get_distance(wfl::WalkerFirstLast)
        return Distance(abs(wfl.w₂.x - wfl.w₁.x),abs(wfl.w₂.y - wfl.w₁.y))
    end

    """
    Takes types Distancae, Variables and Domain and computes the overall boundary conditions dependent distnace of the T6SS run.
    """
    function get_distance(dis::Distance,var::Variables,dom::Domain)
        xboundary,yboundary = get_boundary_distance_functions(var)
        return T6ssDistance(
            xboundary(dis.distance_x,dom.width) + 
            yboundary(dis.distance_y,dom.length))
    end
     
    """
    Takes types Variables, Parameters and Domain and computes the type WalkerFirstLast, then the type Distance and finaly returning type T6ssDistance.
    """
    function get_walker_distance(var::Variables,par::Parameters,dom::Domain)
        return @chain var begin
                    get_walker_position(_,par,dom)
                    get_distance(_)
                    get_distance(_,var,dom)
                end
    end


    """
    Select biological model by specifying the struct: T6ssBiologicalInspired
    Takes types Variables, Parameters and Domain and computes the type WalkerFirstLast, then the type Distance and finaly returning type T6ssDistance.
    """
    function get_walker_distance(::T6ssBiologicalInspired,var::Variables,par::Parameters,dom::Domain)
        return @chain var begin
                    get_walker_position(T6ssBiologicalInspired(),_,par,dom) # Get first and last walker position with one cycle
                    get_distance(_) # get absolute distance for x and y direction
                    get_distance(_,var,dom) # get total distance with consideration of the boundarys
                end
    end




    

    """
    Take types Variables, Parameters, Domain and Iterations and execute the T6SS model such a distribution will emerge dependent on the input parameters.
    """
    function get_experimental_sample_dist_vecs(
        var::Variables,
        par::Parameters,
        dom::Domain,
        iter::Iterations)
        experimental = ThreadsX.map(x -> get_walker_distance(var,par,dom),1:10^iter.simulation_iterations)
        sample = get_sample_distribution(var,dom,iter)
    
        return DistanceExperSampleStruct(
            map(x -> x.distance, experimental),
            sample)
    end



    """
    Select biological model by specifying the struct: T6ssBiologicalInspired
    Take types Variables, Parameters, Domain and Iterations and execute the T6SS model such a distribution will emerge dependent on the input parameters.
    """
    function get_experimental_sample_dist_vecs(
        ::T6ssBiologicalInspired,
        var::Variables,
        par::Parameters,
        dom::Domain,
        iter::Iterations)
        experimental = ThreadsX.map(x -> get_walker_distance(T6ssBiologicalInspired(),var,par,dom),1:10^iter.simulation_iterations)
        sample = get_sample_distribution(var,dom,iter)
    
        return DistanceExperSampleStruct(
            map(x -> x.distance, experimental),
            sample)
    end



    """
    Take type `Domain` and type `Variables` and returns the boundary associated count from 0 to end.
    1. no_boundary
    2. double_periodic
    3. no_flux
    4. cell_boundary
    """
    function  set_range_iter(dom::Domain,var::Variables)
        possibile_options = (
            collect(0:Int((dom.width + dom.length))),
            collect(0:Int((1/2)*(dom.width + dom.length))),
            collect(0:Int((1/2)*(dom.width + (dom.length)))),
            collect(0:Int((1/2)*(dom.width + (dom.length))))
        )
        return boundary_conditional(var.boundary,possibile_options)
    end


    """
    Take type `Domain` and type `Variables` and returns the boundary associated count from 0 to end.
    1. no_boundary
    2. double_periodic
    3. no_flux
    4. cell_boundary
    """
    function  set_range_iter(dom::Domain,var::Variables)
        possibile_options = (
            collect(0:var.Δx:(dom.width + dom.length)),
            collect(0:var.Δx:(1/2)*(dom.width + dom.length)),
            collect(0:var.Δx:((1/2)*(dom.width) + dom.length)),
            collect(0:var.Δx:((1/2)*(dom.width) + dom.length))
        )
        return boundary_conditional(var.boundary,possibile_options)
    end


    







    """
    Takes a vector of reals and returns the count per 0 to a boundary associated end point.
    """
    function get_count(vect_to_count,range)
        counts = [Base.count(i->i==j,vect_to_count) for j in range]
        return counts
    end

