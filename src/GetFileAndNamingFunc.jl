###############################
#### File system functions ####
###############################
//
    
    """
    Setting the struc types for the purpose of saving data to JSON.
    """
    StructTypes.StructType(::Type{SolutionStructure}) = StructTypes.Struct()
    StructTypes.StructType(::Type{DistanceStructure}) = StructTypes.Struct()
    StructTypes.StructType(::Type{DistanceExperSampleStruct}) = StructTypes.Struct()
    StructTypes.StructType(::Type{SolutionVarParDom}) = StructTypes.Struct()
    StructTypes.StructType(::Type{Variables}) = StructTypes.Struct()
    StructTypes.StructType(::Type{Iterations}) = StructTypes.Struct()
    StructTypes.StructType(::Type{Parameters}) = StructTypes.Struct()
    #StructTypes.StructType(::Type{Domain{Int64}}) = StructTypes.Struct()
    StructTypes.StructType(::Type{Domain}) = StructTypes.Struct()

    
    
    
    """
        find(src_dir::String,mindepth::Int64,
        maxdepth::Int64,ext::String;force=false::Bool,
        recursive=false::Bool)

        Find files at multiple depths pased on a pattern.
    """
    function find(src_dir::String,mindepth::Int64,
        maxdepth::Int64,ext::String;force=false::Bool,
        recursive=false::Bool)
        #src_dir::String - Source parent folder for files
        #dst_dir::String - Destination folder for copied renamed files
        #mindepth::Int - Minimum depth for find files
        #maxdepth::Int - Maximum dpeth to find files
        #ext::String - Extension of file to look for 
        dest_fn = "$(src_dir)/list_the_files.txt"
        run(
        pipeline(
        `find $src_dir -mindepth $mindepth -maxdepth $maxdepth -iname "*$ext*"`,
        stdout= dest_fn
        )
        )


        src = @chain dest_fn begin
        open(_)
        readlines(_)
        end

        rm(dest_fn; force=force, recursive=recursive)
        return src
    end



    """
    Take solution of type SolutionVarParDom, which contains the variables, parameters, iterations, domain and solution as input.
    """
    function write_solution_to_file(solution::SolutionVarParDom,path::String)
        # Write solutoin to JSON3
        open(path, "w") do io
            JSON3.pretty(io, solution,allow_inf=true)
        end
    end

    

    """
    Take `String` representing the data file path, along with the `DataType` struct and loads data to inputted `DataType`.
    """
    function read_solution_from_memory(path::String,json_struct::DataType) 
        # Read Jason string 
        json_string = read(path, String)
        # Convert to data format
        input = JSON3.read(json_string,json_struct)
        return input
    end



    """
    Takes an `Int64` and `String` as input and returns a `String` in the form "Figure*`Int`*`String`". 
    As an example, say `join_figure_number_letter(1,"a")` will return "Figure1a". 
    Acts as part of the label for saving figures and data.
    """
    join_figure_number_letter(fig_numeric_label::Int64,fig_alpha_label::String) = "Figure"*string(fig_numeric_label)*fig_alpha_label
    
    """
    Takes a `Tuple` or a `String` and returns the `join` function. Used when a dataset covers multiple figures. Example ("a","b","c") = "abc"
    """
    join_name_letters(fig_alpha_label::Union{Tuple,String}) = join(fig_alpha_label)

    """
        Takes the `Int64` related to the figure number and returns a path type string. This function is used in saving figures to a specific folder dependent on the figure number.
    """
    append_forward_slash_figure(fig_numeric_label::Int64) = "Figure"*string(fig_numeric_label)*"/"

    """
        Takes `Int64` numeric label for the figure and either a `Tuple` of `String` or just `String` and returns as example "Figure1abc". 
    """
    join_figure_number_tuple_letters(fig_numeric_label::Int64,fig_alpha_label::Union{Tuple,String}) = "Figure"*string(fig_numeric_label)*join_name_letters(fig_alpha_label)

    """
        Takes the path to be used for the data, called `data_root`, the `Int64` numeric label and either a `Tuple` of `String` or just `String` to return the full path for the data to be saved  
    """
    data_path(data_root::String,fig_numeric_label::Int64,fig_alpha_label::Union{Tuple,String}) = data_root*"SolutionData"*join_figure_number_tuple_letters(fig_numeric_label,fig_alpha_label)*".json"

    
//


