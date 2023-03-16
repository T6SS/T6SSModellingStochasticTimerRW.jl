
    function generate_t6ss_data(input_values::T6ssDataGenerateInput)
        @unpack Δt, Δx, T, X, aspect,boundary,simulation_iterations,parameter_iterations,λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root = 
            input_values
        var = Variables(Δt, Δx, T, X, aspect,boundary)
        par =  Parameters(λ₀,λ₁,h)
        dom = dimensions(var)
        iter = Iterations(simulation_iterations,parameter_iterations)
        λ₀_vec = range(0,λmax,iter.parameter_iterations)
        h_vec = range(0,hmax,iter.parameter_iterations)
        a_vec = range(0,amax,iter.parameter_iterations)
        D_vec = range(dif_min,dif_max,iter.parameter_iterations) 

        
        

        # Set data path and file names
        fig_numeric_label =(3,4,5,6)
        parameter_span = "Span"
        
        data_path_json_vec = [
            data_root*
            join_figure_number_letter(x,parameter_span)*
            "ParamIter"*
            string(z) 
                for x in fig_numeric_label 
                for z in 1:iter.parameter_iterations] .* 
                ".json"

                
        

        # Compute solutions of parameters
        # λ₀_vec
        parλ₀ = map(λ₀_vec) do io
            @set par.λ₀ = io
        end
        
        solutionλ₀ = map(parλ₀) do p
            Δp = Δ(p,var)       
            res = get_experimental_sample_dist_vecs(var, Δp, dom, iter)
            solu = SolutionVarParDom(var,p,dom,iter,res.experimental,res.sample)
        end
        
        # h_vec
        parh = map(h_vec) do io
            @set par.h = io
        end
        
        solutionh = map(parh) do p
            Δp = Δ(p,var)       
            res = get_experimental_sample_dist_vecs(var, Δp, dom, iter)
            solu = SolutionVarParDom(var,p,dom,iter,res.experimental,res.sample)
        end

        #a_vec
        vara = map(a_vec) do io
            @set var.aspect = io
        end
        
        solutiona = map(vara) do v       
            dom = dimensions(v)
            Δp = Δ(par,v)       
            res = get_experimental_sample_dist_vecs(v, Δp, dom, iter)
            solu = SolutionVarParDom(v,par,dom,iter,res.experimental,res.sample)
        end
        
        # D_vec
        Δx = 0.1
        aspect = 1.6
        X = 2
        var = @set var.aspect = aspect
        dom = dimensions(var)
        parD = map(D_vec) do io    
            @set par.h = io
        end

        
        solutionD = map(parD) do io
            parDLitRate = @set io.h = round(LitRate(io.h,Δx),digits=0)
            Δp = Δ(parDLitRate,var)       
            res = get_experimental_sample_dist_vecs(var, Δp, dom, iter) * Δx 
        end

        var = @set var.Δx = Δx
        var = @set var.X = X
        dom = dimensions(var)
        sols = [SolutionVarParDom(var,parD[i],dom,iter,solutionD[i].experimental,solutionD[i].sample) for i in 1:iter.parameter_iterations]


        # Append and write solutions
        solution_params = vcat(solutionλ₀,solutionh,solutiona,sols)
        write_solution_to_file.(solution_params,data_path_json_vec)
    end


