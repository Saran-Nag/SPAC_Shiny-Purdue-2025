from shiny import ui, render, reactive
import anndata as ad
import pandas as pd
import spac.visualization


def scatterplot_server(input, output, session, shared):
    @reactive.Calc
    def get_scatterplot_names():
        if (shared['obsm_names'].get() is not None and 
            shared['var_names'].get() is not None):
            obsm_list = shared['obsm_names'].get()
            var_list = shared['var_names'].get()
            obsm_dict = {item: item for item in obsm_list}
            features_dict = {item: item for item in var_list}
            dict = {"Annotated Tables" : obsm_dict, "Features" : features_dict}

            return dict
        return []

    @reactive.Effect
    def update_select_input_layer_scatter():
        choices = get_scatterplot_names()
        ui.update_select("scatter_x", choices=choices)
        ui.update_select("scatter_y", choices=choices)
        return

    @reactive.Calc
    def get_scatterplot_coordinates_x():
        adata = ad.AnnData(
            X=shared['X_data'].get(), 
            var=pd.DataFrame(shared['var_data'].get()), 
            obsm=shared['obsm_data'].get(), 
            layers=shared['layers_data'].get()
        )
        obsm = shared['obsm_names'].get()
        features = shared['var_names'].get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_x()

        if selection in obsm:
            coords = adata.obsm[selection]

            # Extract the first column for x-coordinates
            x_coords = coords[:, 0]  
            return x_coords
        elif selection in features and layer_selection == "Original":
            column_index = adata.var_names.get_loc(selection)

            # Extract the column corresponding to the feature
            x_coords = adata.X[:, column_index]  
            return x_coords
        elif selection in features and layer_selection != "Original":
            column_index = adata.var_names.get_loc(selection)
            new_layer = adata.layers[layer_selection]

            # Extract the column corresponding to the feature
            x_coords = new_layer[:, column_index]
            return x_coords

        return None



    @reactive.Calc
    def get_scatterplot_coordinates_y():
        adata = shared['adata_main'].get()
        obsm = shared['obsm_names'].get()
        features = shared['var_names'].get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_y()

        if selection in obsm:
            coords = adata.obsm[selection]

            # Extract the second column for y-coordinates
            y_coords = coords[:, 1]  
            return y_coords
        elif selection in features and layer_selection == "Original":
            column_index = adata.var_names.get_loc(selection)

            # Extract the column corresponding to the feature
            y_coords = adata.X[:, column_index]  
            return y_coords
        elif selection in features and layer_selection != "Original":
            column_index = adata.var_names.get_loc(selection)
            new_layer = adata.layers[layer_selection]

            # Extract the column corresponding to the feature
            y_coords = new_layer[:, column_index]  
            return y_coords

        return None


    # Track the UI state for scatterplot dropdowns
    scatter_ui_initialized = reactive.Value(False)

    @reactive.effect
    def scatter_reactivity():
        btn = input.scatter_color_check()
        if btn and not scatter_ui_initialized.get():
            # Insert the color selection dropdown if not already initialized
            dropdown = ui.input_select(
                "scatter_color", 
                "Select Feature", 
                choices=shared['var_names'].get()
            )
            ui.insert_ui(
                ui.div({"id": "inserted-scatter_dropdown"}, dropdown),
                selector="#main-scatter_dropdown",
                where="beforeEnd",
            )
            scatter_ui_initialized.set(True)
        elif not btn and scatter_ui_initialized.get():
            # Remove the color selection dropdown if it exists
            ui.remove_ui("#inserted-scatter_dropdown")
            scatter_ui_initialized.set(False)

    @reactive.Calc
    def get_color_values():
        selected_feature = input.scatter_color()
        if selected_feature is None:
            return None
        adata = ad.AnnData(
            X=shared['X_data'].get(), 
            var=pd.DataFrame(shared['var_data'].get())
        )
        if selected_feature in adata.var_names:
            column_index = adata.var_names.get_loc(selected_feature)
            color_values = adata.X[:, column_index]
            return color_values
        return None 

    @output
    @render.plot
    @reactive.event(input.go_scatter, ignore_none=True)
    def spac_Scatter():
        x_points = get_scatterplot_coordinates_x()
        y_points = get_scatterplot_coordinates_y()
        btn = input.scatter_color_check()
        x_label = input.scatter_x()
        y_label = input.scatter_y()
        title = f"Scatterplot: {x_label} vs {y_label}"
        if btn is False:
            fig, ax = spac.visualization.visualize_2D_scatter(x_points,y_points)
            ax.set_title(title, fontsize=14)
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            return ax
        elif btn is True:
            fig1, ax1 = spac.visualization.visualize_2D_scatter(
                x_points,y_points, 
                labels=get_color_values()
            )
            ax1.set_title(title, fontsize=14)
            ax1.set_xlabel(x_label)
            ax1.set_ylabel(y_label)
            for color in fig1.axes:
                if hasattr(color, "get_ylabel") and color != ax1:
                    color.set_ylabel(f"Colored by: {input.scatter_color()}")
            return ax1
