#####################################################
#### Visualise the multiple iterations solutions ####
#####################################################
//
    # define data in terms of a vector of random Walkers
    # for use in algebra of graphics
    AlgebraOfGraphics.data(w::Vector{Walker2D{Int64}}) = RandomWalker.position.(w)




    # take solution of T6SS
    # visualise solution
    # copy parameters used to figure folder for reference
    function visualise_save_solution(
            solution::NamedTuple,
            fig_num::Int,
            fig_sub_num::Union{Tuple,String},
            fig_folder::String,
            parameter_file::String)
        figure_names = fig_labels.(fig_num,fig_sub_num)
        figure_folder = fig_labels(fig_num)
        data_name = fig_labels(fig_num,fig_sub_num)
        
        # Get maximum distance size based on 
        # xwidth and ylength of domain
        sindx = size(solution[:indx],1)
        
        # Get a dark cyan color for bar plot
        dark_cyan = colorant"#006390"

        # Set axis sizes
        lims = Tuple(reshape([-1/2,1/2].*[xwidth ylength],4,1))
        #axishtmp = (width = 512, height = 512,aspect=1,limits = lims)
        axis = (width = 512, height = 512,aspect=1,xlimits = (0,sindx))
        
        # Set figure resolution
        #size_inches = (4, 3)
        #size_pt = 72 .* size_inches
        
        # Plot heatmap
        #=
        htmp = 
            data(solution[:machine_positions]) * 
            mapping("1" => L"x","2" => L"y") * 
            AlgebraOfGraphics.histogram()
        fhtmp = draw(htmp;axis = axishtmp,colorbar=(position=:top, size=25, label="Count"))
        =#
        

        # Plot distances for sample and experimental
        tab_dist_count = (
            x = repeat(solution[:indx],outer = sindx),
            y = vcat(
                solution[:count_experimental],
                solution[:count_sample]),
            c = repeat(["Experimental", "Sample"],inner = sindx))
        dist_count_bar = data(tab_dist_count) *
            mapping(:x => L"d",:y => L"count",color = :c => "Group") *
            visual(BarPlot,alpha = 0.4)
        fdist_count_bar = draw(
            dist_count_bar;
            axis = axis, 
            legend=(position=:top, titleposition=:left, framevisible=true, padding=5))
        

        # Plot count norm ratio
        #=
        norm_count_bar = 
            data(solution) *
            mapping(:indx => L"d",:count_norm => L"count") * 
            visual(BarPlot,color=dark_cyan)
        fnorm_count_bar = draw(norm_count_bar;axis = axis)
        =#
        # Store figure in tuple 

        #figures = (fdist_count_bar,fnorm_count_bar)
        figures = fdist_count_bar

        
        

        # Save all figures
        #[save(fig_folder*figure_folder*figure_names[i]*".pdf",figures[i],pt_per_unit = 1) for i in 1:size(figures,1)]
        save(fig_folder*figure_folder*figure_names*".pdf",figures,pt_per_unit = 1)
        # Copy parametes used in simultions to figure folder
        cp(
            parameter_file, 
            fig_folder*figure_folder*"ParameterConfig"*data_name*".txt"; 
            force = true)

    end

    function visualise_save_solution(
            solution::NamedTuple,
            fig_num::Int,
            fig_sub_num::Union{Tuple,String},
            fig_folder::String,
            parameter_file::String;
            walker_ic::String = "random")
        figure_names = fig_labels.(fig_num,fig_sub_num)
        figure_folder = fig_labels(fig_num)
        data_name = fig_labels(fig_num,fig_sub_num)
        
        # Get maximum distance size based on 
        # xwidth and ylength of domain
        sindx = size(solution[:indx],1)
        
        # Get a dark cyan color for bar plot
        dark_cyan = colorant"#006390"

        # Set axis sizes
        lims = Tuple(reshape([-1/2,1/2].*[xwidth ylength],4,1))
        axishtmp = (width = 512, height = 512,aspect=1,limits = lims)
        axis = (width = 512, height = 512,aspect=1,xlimits = (0,sindx))
        
        # Set figure resolution
        #size_inches = (4, 3)
        #size_pt = 72 .* size_inches
        if walker_ic == "random"
        end
        # Plot heatmap
        #=
        if walker_ic == "random"
        htmp = 
            data(solution[:machine_final]) * 
            mapping("1" => L"x","2" => L"y") * 
            AlgebraOfGraphics.histogram()
        else
            htmp = 
            data(solution[:machine_positions]) * 
            mapping("1" => L"x","2" => L"y") * 
            AlgebraOfGraphics.histogram()
        end
        fhtmp = draw(htmp;axis = axishtmp,colorbar=(position=:top, size=25, label="Count"))
        =#

        

        # Plot distances for sample and experimental
        tab_dist_count = (
            x = repeat(solution[:indx],outer = sindx),
            y = vcat(
                solution[:count_experimental],
                solution[:count_sample]),
            c = repeat(["Experimental", "Sample"],inner = sindx))
        dist_count_bar = data(tab_dist_count) *
            mapping(:x => L"d",:y => L"count",color = :c => "Group") *
            visual(BarPlot)
        fdist_count_bar = draw(
            dist_count_bar;
            axis = axis, 
            legend=(position=:top, titleposition=:left, framevisible=true, padding=5))
        
        #=
        # Plot count norm ratio
        norm_count_bar = 
            data(solution) *
            mapping(:indx => L"d",:count_norm => L"count") * 
            visual(BarPlot,color=dark_cyan)
        fnorm_count_bar = draw(norm_count_bar;axis = axis)
        
        # Store figure in tuple 
        =#
        #figures = (fhtmp,fdist_count_bar,fnorm_count_bar)
        figures = fdist_count_bar
        
        

        # Save all figures
        #[save(fig_folder*figure_folder*figure_names[i]*".pdf",figures[i],pt_per_unit = 1) for i in 1:3]
        save(fig_folder*figure_folder*figure_names*".pdf",figures,pt_per_unit = 1)
        # Copy parametes used in simultions to figure folder
        cp(
            parameter_file, 
            fig_folder*figure_folder*"ParameterConfig"*data_name*".txt"; 
            force = true)

    end
//


