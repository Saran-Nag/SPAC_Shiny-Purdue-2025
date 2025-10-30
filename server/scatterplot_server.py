from shiny import ui, render, reactive
import anndata as ad
import pandas as pd
import spac.visualization
import seaborn as sns
from utils.datashader_utils import scatter_heatmap
import matplotlib.pyplot as plt


def scatterplot_server(input, output, session, shared):
    @reactive.calc
    def get_scatterplot_names():
        obsm_list = shared['obsm_names'].get()
        var_list = shared['var_names'].get()

        if obsm_list is None or var_list is None:
            return []

        return {
            "Annotated Tables": {item: item for item in obsm_list},
            "Features": {item: item for item in var_list},
        }

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

        obsm_names = shared['obsm_names'].get()
        features = shared['var_names'].get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_x()

        if selection in obsm_names:
            return adata.obsm[selection][:, 0]

        if selection in features:
            col_idx = adata.var_names.get_loc(selection)
            if layer_selection == "Original":
                return adata.X[:, col_idx]
            else:
                return adata.layers[layer_selection][:, col_idx]

        return None



    @reactive.Calc
    def get_scatterplot_coordinates_y():
        adata = shared['adata_main'].get()
        obsm_names = shared['obsm_names'].get()
        features = shared['var_names'].get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_y()

        if selection in obsm_names:
            return adata.obsm[selection][:, 1]

        if selection in features:
            col_idx = adata.var_names.get_loc(selection)
            if layer_selection == "Original":
                return adata.X[:, col_idx]
            else:
                return adata.layers[layer_selection][:, col_idx]

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
        x = get_scatterplot_coordinates_x()
        y = get_scatterplot_coordinates_y()
        color_enabled = input.scatter_color_check()
        heatmap_mode = input.scatter_heatmap_mode()
        x_label = input.scatter_x()
        y_label = input.scatter_y()
        title = f"Scatterplot: {x_label} vs {y_label}"
        font_size = input.scatter_font_size()
        with sns.plotting_context(rc={"font.size": font_size,
                                      "axes.labelsize": font_size,
                                      "xtick.labelsize": font_size,
                                      "ytick.labelsize": font_size,
                                      "legend.fontsize": font_size,
                                      "axes.titlesize": font_size * 1.2
        }):

        if heatmap_mode:
            color = get_color_values() if color_enabled else None
            img = scatter_heatmap(x, y, color)
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.imshow(img, aspect='auto')
            ax.set_title(title, fontsize=14)
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.axis('on')  # Show axes
            return fig
        else:
            if color_enabled:
                fig, ax = spac.visualization.visualize_2D_scatter(
                    x, y, labels=get_color_values()
                )
                for a in fig.axes:
                    if hasattr(a, "get_ylabel") and a != ax:
                        a.set_ylabel(f"Colored by: {input.scatter_color()}")
            else:
                fig, ax = spac.visualization.visualize_2D_scatter(x, y)

            ax.set_title(title, fontsize=font_size * 1.2)
            ax.set_xlabel(x_label, fontsize=font_size)
            ax.set_ylabel(y_label, fontsize=font_size)

            return ax
