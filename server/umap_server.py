from shiny import ui, render, reactive
import anndata as ad
import pandas as pd
import spac.visualization


def umap_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_umap1, ignore_none=True)
    def spac_UMAP():
        adata = ad.AnnData(
            X=shared['X_data'].get(), 
            var=pd.DataFrame(shared['var_data'].get()), 
            obsm=shared['obsm_data'].get(), 
            obs=shared['obs_data'].get(), 
            dtype=shared['X_data'].get().dtype, 
            layers=shared['layers_data'].get()
        )
        point_size=input.umap_slider_1()
        if adata is not None:
            if input.umap_rb() == "Feature":
                if input.umap_layer() == "Original":
                    layer = None
                else:
                    layer = input.umap_layer()
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata, 
                    method=input.plottype(), 
                    feature=input.umap_rb_feat(), 
                    layer=layer, 
                    point_size=point_size
                )
                ax.set_title(
                    f"{input.plottype().upper()}: {input.umap_rb_feat()}",
                    fontsize=14
                )
                ax.set_xlabel(f"{input.plottype().upper()} 1")
                ax.set_ylabel(f"{input.plottype().upper()} 2")
                for color_ax in fig.axes:
                    if hasattr(color_ax, "get_ylabel") and color_ax != ax:
                        color_ax.set_ylabel(
                            f"Colored by: {input.umap_rb_feat().upper()}",
                            fontsize=12
                        )
                return fig
            elif input.umap_rb() == "Annotation":
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata, 
                    method=input.plottype(), 
                    annotation=input.umap_rb_anno(), 
                    point_size=point_size
                )
                ax.set_title(
                    f"{input.plottype().upper()}: {input.umap_rb_anno()}",
                    fontsize=14
                )
                ax.set_xlabel(f"{input.plottype().upper()} 1")
                ax.set_ylabel(f"{input.plottype().upper()} 2")
                return fig
        return None

    # Track the UI state
    umap_annotation_initialized = reactive.Value(False)
    umap_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.umap_rb()

            if btn == "Annotation":
                if not umap_annotation_initialized.get():
                    # Create the Annotation dropdown
                    dropdown = ui.input_select(
                        "umap_rb_anno", 
                        "Select an Annotation", 
                        choices=shared['obs_names'].get(),
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_anno"}, dropdown),
                        selector="#main-ump_rb_dropdown_anno",
                        where="beforeEnd",
                    )
                    # Update the state
                    umap_annotation_initialized.set(True)
                # Remove the Feature dropdown and table
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)

            elif btn == "Feature":
                if not umap_feature_initialized.get():
                    # Create the Feature dropdown
                    dropdown1 = ui.input_select(
                        "umap_rb_feat", 
                        "Select a Feature", 
                        choices=shared['var_names'].get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat",
                        where="beforeEnd",
                    )

                    # Create the Table dropdown
                    new_choices = shared['layers_names'].get() + ["Original"]
                    table_umap = ui.input_select(
                        "umap_layer", 
                        "Select a Table", 
                        choices=new_choices, 
                        selected=["Original"]
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-umap_table"}, table_umap),
                        selector="#main-ump_table_dropdown_feat",
                        where="beforeEnd",
                    )
                    # Update the state
                    umap_feature_initialized.set(True)
                # Remove the Annotation dropdown
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)

            elif btn == "None":
                # Remove all dropdowns and reset states
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)


    @output
    @render.plot
    @reactive.event(input.go_umap2, ignore_none=True)
    def spac_UMAP2():
        adata = ad.AnnData(
            X=shared['X_data'].get(), 
            var=pd.DataFrame(shared['var_data'].get()), 
            obsm=shared['obsm_data'].get(), 
            obs=shared['obs_data'].get(), 
            dtype=shared['X_data'].get().dtype, 
            layers=shared['layers_data'].get()
        )
        point_size_2=input.umap_slider_2()
        if adata is not None:
            if input.umap_rb2() == "Feature":
                if input.umap_layer2() == "Original":
                    layer2 = None
                else:
                    layer2 = input.umap_layer2()
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata, 
                    method=input.plottype2(), 
                    feature=input.umap_rb_feat2(), 
                    layer=layer2, 
                    point_size=point_size_2
                )
                ax.set_title(
                    f"{input.plottype2().upper()}: {input.umap_rb_feat2()}",
                    fontsize=14
                )
                ax.set_xlabel(f"{input.plottype2().upper()} 1")
                ax.set_ylabel(f"{input.plottype2().upper()} 2")
                for color_ax in fig.axes:
                    if hasattr(color_ax, "get_ylabel") and color_ax != ax:
                        color_ax.set_ylabel(
                            f"Colored by: {input.umap_rb_feat2()}", 
                            fontsize=12
                        )
                return fig
            elif input.umap_rb2() == "Annotation":
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata, 
                    method=input.plottype2(), 
                    annotation=input.umap_rb_anno2(), 
                    point_size=point_size_2
                )
                ax.set_title(
                    f"{input.plottype2().upper()}: {input.umap_rb_anno2()}", fontsize=14
                )
                ax.set_xlabel(f"{input.plottype2().upper()} 1")
                ax.set_ylabel(f"{input.plottype2().upper()} 2")
                return fig
        return None

    # Track the UI state
    umap2_annotation_initialized = reactive.Value(False)
    umap2_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity2():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.umap_rb2()

            if btn == "Annotation":
                if not umap2_annotation_initialized.get():
                    dropdown = ui.input_select(
                        "umap_rb_anno2", 
                        "Select an Annotation", 
                        choices=shared['obs_names'].get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_anno2"}, dropdown),
                        selector="#main-ump_rb_dropdown_anno2",
                        where="beforeEnd",
                    )
                    umap2_annotation_initialized.set(True)
                if umap2_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat2")
                    ui.remove_ui("#inserted-umap_table2")
                    umap2_feature_initialized.set(False)

            elif btn == "Feature":
                if not umap2_feature_initialized.get():
                    dropdown1 = ui.input_select(
                        "umap_rb_feat2", 
                        "Select a Feature", 
                        choices=shared['var_names'].get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat2"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat2",
                        where="beforeEnd",
                    )

                    new_choices = shared['layers_names'].get() + ["Original"]
                    table_umap_1 = ui.input_select(
                        "umap_layer2", "Select a Table", choices=new_choices, selected=["Original"]
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-umap_table2"}, table_umap_1),
                        selector="#main-ump_table_dropdown_feat2",
                        where="beforeEnd",
                    )
                    umap2_feature_initialized.set(True)
                if umap2_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno2")
                    umap2_annotation_initialized.set(False)

            elif btn == "None":
                if umap2_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno2")
                    umap2_annotation_initialized.set(False)
                if umap2_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat2")
                    ui.remove_ui("#inserted-umap_table2")
                    umap2_feature_initialized.set(False)
