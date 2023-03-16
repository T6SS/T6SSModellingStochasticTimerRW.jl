
    function solve_iterations()
        include("parameters/ParameterConfigFigure3d.jl")
        # λ₀Δt
        d = map(it -> mean(
                    get_distance(
                        get_walker_position(it,λ₁Δt,hΔt,t⃗,N,xwidth,ylength,boundary,boundary_conditional),
                        boundary)),λ₀Δt)
        solution = (λ₀Δt = λ₀Δt, mean_distances = d)
        # Create json path
        json_path = data_path(data_root,fig_num,fig_sub_num)    
        # Write to file
        write_solution_to_file(solution,json_path)
        
        figure_names = fig_labels.(fig_num,fig_sub_num)
        figure_folder = fig_labels(fig_num)
        data_name = fig_labels(fig_num,fig_sub_num)
        dark_cyan = colorant"#006390"
        axis = (width = 512, height = 512,aspect=1)
        x = solution[1]
        y = solution[2]
        df = (; x, y)
        specs = data(df) * 
        mapping(:x => L"\lambda_0", :y => L"d") * 
        visual(Scatter,color=dark_cyan)
        figure = draw(specs,axis = axis)
        save(fig_folder*figure_folder*figure_names*".pdf",figure,pt_per_unit = 1)

        cp(
            "parameters/ParameterConfigFigure3d.jl", 
            fig_folder*figure_folder*"ParameterConfig"*data_name*".txt"; 
            force = true)





        include("parameters/ParameterConfigFigure4d.jl")
        # hΔt
        d = map(it -> mean(
                    get_distance(
                        get_walker_position(λ₀Δt,λ₁Δt,it,t⃗,N,xwidth,ylength,boundary,boundary_conditional),
                        boundary)),hΔt)
        solution = (hΔt = hΔt, mean_distances = d)
        # Create json path
        json_path = data_path(data_root,fig_num,fig_sub_num)    
        # Write to file
        write_solution_to_file(solution,json_path)

        
        figure_names = fig_labels.(fig_num,fig_sub_num)
        figure_folder = fig_labels(fig_num)
        data_name = fig_labels(fig_num,fig_sub_num)
        dark_cyan = colorant"#006390"
        axis = (width = 512, height = 512,aspect=1)
        x = solution[1]
        y = solution[2]
        df = (; x, y)
        specs = data(df) * 
        mapping(:x => L"h", :y => L"d") * 
        visual(Scatter,color=dark_cyan)
        figure = draw(specs,axis = axis)
        save(fig_folder*figure_folder*figure_names*".pdf",figure,pt_per_unit = 1)

        cp(
            "parameters/ParameterConfigFigure4d.jl", 
            fig_folder*figure_folder*"ParameterConfig"*data_name*".txt"; 
            force = true)

    
    
    
    
        include("parameters/ParameterConfigFigure5d.jl")                
        # ylength
        d = map(it -> mean(
                    get_distance(
                        get_walker_position(λ₀Δt,λ₁Δt,hΔt,t⃗,N,xwidth,it,boundary,boundary_conditional),
                        boundary)),ylength)
        solution = (ylength = ylength, mean_distances = d)
        # Create json path
        json_path = data_path(data_root,fig_num,fig_sub_num)    
        # Write to file
        write_solution_to_file(solution,json_path)
    
        
        figure_names = fig_labels.(fig_num,fig_sub_num)
        figure_folder = fig_labels(fig_num)
        data_name = fig_labels(fig_num,fig_sub_num)
        dark_cyan = colorant"#006390"
        axis = (width = 512, height = 512,aspect=1)
        x = solution[1]
        y = solution[2]
        df = (; x, y)
        specs = data(df) * 
        mapping(:x => L"length", :y => L"d") * 
        visual(Scatter,color=dark_cyan)
        figure = draw(specs,axis = axis)
        save(fig_folder*figure_folder*figure_names*".pdf",figure,pt_per_unit = 1)

        cp(
            "parameters/ParameterConfigFigure5d.jl", 
            fig_folder*figure_folder*"ParameterConfig"*data_name*".txt"; 
            force = true)

    
    
    
        include("parameters/ParameterConfigFigure6d.jl")                
    
        # ylength
        d = map(it -> mean(
                    get_distance(
                        get_walker_position(λ₀Δt,λ₁Δt,it,t⃗,N,xwidth,ylength,boundary,boundary_conditional),
                        boundary)),hΔt)
        solution = (DΔt = hΔt.*10, mean_distances = d)
        # Create json path
        json_path = data_path(data_root,fig_num,fig_sub_num)    
        # Write to file
        write_solution_to_file(solution,json_path)


        figure_names = fig_labels.(fig_num,fig_sub_num)
        figure_folder = fig_labels(fig_num)
        data_name = fig_labels(fig_num,fig_sub_num)
        dark_cyan = colorant"#006390"
        axis = (width = 512, height = 512,aspect=1)
        x = solution[1]
        y = solution[2]
        df = (; x, y)
        specs = data(df) * 
        mapping(:x => L"D", :y => L"d") * 
        visual(Scatter,color=dark_cyan)
        figure = draw(specs,axis = axis)
        save(fig_folder*figure_folder*figure_names*".pdf",figure,pt_per_unit = 1)

        cp(
            "parameters/ParameterConfigFigure6d.jl", 
            fig_folder*figure_folder*"ParameterConfig"*data_name*".txt"; 
            force = true)



    end
    
//






