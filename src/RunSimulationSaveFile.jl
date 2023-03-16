############################################
#### Solve, visualise and save solution ####
############################################
//
    # aggregate all functions into this functions
    # all variables/parameters are defined in the 
    # parameter config file
    function run_simulation(parameter_file::String)
        # Solve 
        solution = solve(parameter_file)


        # Store visuals
        visualise_save_solution(
            solution,fig_num,fig_sub_num,fig_folder,parameter_file)
            
        # Create json path
        json_path = data_path(data_root,fig_num,fig_sub_num)    
        # Write to file
        write_solution_to_file(solution,json_path)
        return json_path
    end

    function run_simulation(parameter_file::String,walker_ic::String = "random")
        # Solve 
        if walker_ic == "random"
            solution = solve(parameter_file; walker_ic = walker_ic)
        else
            solution = solve(parameter_file)
        end


        # Store visuals
        visualise_save_solution(
            solution,fig_num,fig_sub_num,fig_folder,parameter_file)
            
            # Create json path
        json_path = data_path(data_root,fig_num,fig_sub_num)    
        # Write to file
        write_solution_to_file(solution,json_path)
        return json_path
    end
//
