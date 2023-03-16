"""
    rand_walker(dom::Domain,var::Variables)

    Get a walker position randomly in the domain.
    Use Δx as the step.
"""
    function rand_walker(dom::Domain,var::Variables)
        return Walker2D(
            Int(rand(range(dom.width,var))),
            Int(rand(range(dom.length,var)))
            )
    end
   
    

    function boundary_conditional(boundary::String,possibile_options::Tuple)
        if boundary == "no_boundary"
            possibile_options[1]
        elseif boundary == "double_periodic"
            possibile_options[2]
        elseif boundary == "no_flux"
            possibile_options[3]
        elseif boundary == "cell_boundary"
            possibile_options[4]
        else
            error("boundary has to be no_boundary,double_periodic,no_flux or cell_boundary, not $(boundary)")
        end
    end
    
    function get_boundary_distance_functions(var::Variables)
        if var.boundary == "no_boundary"
            return (xboundary = get_no_periodic_boundary_distance, yboundary = get_no_periodic_boundary_distance)
        elseif var.boundary == "double_periodic"
            return (xboundary = get_periodic_boundary_distance, yboundary = get_periodic_boundary_distance)
        elseif var.boundary == "no_flux"
            return (xboundary = get_periodic_boundary_distance, yboundary = get_no_periodic_boundary_distance)
        elseif var.boundary == "cell_boundary"
            return (xboundary = get_periodic_boundary_distance, yboundary = get_no_periodic_boundary_distance)
        else
            error("boundary = $(boundary) is not no_boundary, double_periodic, no_flux or cell_boundary.")
        end
    end
    

    boundary_conditions(var::Variables) = boundary_conditional(
                    var.boundary,(
                        RandomWalker.update,
                        RandomWalker.updateperiod,
                        RandomWalker.updatenoflux,
                        RandomWalker.updatecell))
    
    Telegraph.update(state,par::Parameters) = Telegraph.update(state,par.λ₀,par.λ₁)
    
    function walker_update(
        walker::Walker2D,
        state::Int64,
        var::Variables,
        par::Parameters,
        dom::Domain
        )
        rw_update = boundary_conditions(var)
        return rw_update(walker, state,par.h, dom.width, dom.length)        
    end
    

    """
        T6ss_next_choice(par::Parameters,r_uni::Float64)
        Get the next step or state from a conditional
    """
    function T6ss_next_choice(par::Parameters,r_uni::Float64)
        if r_uni <= par.λ₁
            # Assumed that state = 1 and moving to state = 0
            return state = 0
        elseif r_uni <= par.h + par.λ₁
            return Step2D(0,-1)
        elseif  r_uni <= 2*par.h + par.λ₁
            return Step2D(0,1)
        elseif  r_uni <= 3*par.h + par.λ₁
            return Step2D(-1,0)
        elseif  r_uni <= 4*par.h + par.λ₁
            return Step2D(1,0)
        elseif  r_uni <= (4*par.h + par.λ₁+ par.λ₀)
            return state = 1
        elseif r_uni <= (1 - (4*par.h + par.λ₀ + par.λ₁))
            return Step2D(0,0)
        else
            # Nothig left
        end
    end