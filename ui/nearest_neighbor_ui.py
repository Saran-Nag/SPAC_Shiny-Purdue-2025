"""
Nearest neighbor visualization UI module for SPAC Shiny application.

This module provides the user interface components for visualizing precomputed
nearest neighbor distances using the visualize_nearest_neighbor functionality.
"""

from shiny import ui


def nearest_neighbor_ui():
    """
    Create the nearest neighbor visualization UI.
    
    Returns
    -------
    shiny.ui.NavPanel
        UI components for the nearest neighbor visualization feature
    """
    return ui.nav_panel(
        "Nearest Neighbors",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        3,
                        ui.div(
                            ui.h4("Core Parameters",
                                  class_="accessible-heading"),
                            
                            # Core functionality parameters
                            ui.input_select(
                                "nn_source_label",
                                "Source Anchor Cell Label",
                                choices=[]
                            ),
                            ui.input_selectize(
                                "nn_target_label",
                                "Target Cell Label",
                                choices=[],
                                multiple=True,
                                options={
                                    "placeholder":
                                    "Select targets (empty for 'All')"
                                }
                            ),
                            ui.input_select(
                                "nn_image_id",
                                "ImageID",
                                choices=["None"]
                            ),
                            
                            ui.hr(),
                            
                            # Plot configuration
                            ui.h5("Plot Configuration",
                                  class_="accessible-heading"),
                            ui.input_select(
                                "nn_plot_method",
                                "Plot Method",
                                choices=["numeric", "distribution"],
                                selected="numeric"
                            ),
                            ui.panel_conditional(
                                "input.nn_plot_method === 'numeric'",
                                ui.input_select(
                                    "nn_plot_type_numeric",
                                    "Plot Type",
                                    choices=["boxen", "box", "violin",
                                             "strip", "swarm"],
                                    selected="boxen"
                                ),
                            ),
                            ui.panel_conditional(
                                "input.nn_plot_method === 'distribution'",
                                ui.input_select(
                                    "nn_plot_type_distribution",
                                    "Plot Type",
                                    choices=["kde", "hist", "ecdf"],
                                    selected="kde"
                                ),
                            ),
                            ui.input_checkbox(
                                "nn_log_scale",
                                "Log Scale",
                                value=False
                            ),
                            ui.input_checkbox(
                                "nn_facet_plot",
                                "Facet Plot",
                                value=False
                            ),
                        ),
                        
                        ui.div(
                            ui.h5("Advanced Settings",
                                  class_="accessible-heading"),
                            
                            ui.input_text(
                                "nn_distance_table",
                                "Nearest Neighbor Associated Table",
                                value="spatial_distance"
                            ),
                            ui.input_text(
                                "nn_color_mapping",
                                "Defined Color Mapping",
                                value="None",
                                placeholder="Enter color map name or 'None'"
                            ),
                            
                            ui.hr(),
                            
                            # Figure configuration parameters
                            ui.h5("Figure Configuration",
                                  class_="accessible-heading"),
                            ui.input_numeric(
                                "nn_x_axis_rotation",
                                "X Axis Label Rotation (degrees)",
                                value=0,
                                min=-90,
                                max=90
                            ),
                            ui.input_checkbox(
                                "nn_shared_x_title",
                                "Shared X Axis Title",
                                value=True
                            ),
                            ui.input_numeric(
                                "nn_x_title_fontsize",
                                "X Axis Title Font Size",
                                value=12,
                                min=8,
                                max=20
                            ),
                            
                            ui.input_numeric(
                                "nn_figure_width",
                                "Figure Width",
                                value=12,
                                min=4,
                                max=20
                            ),
                            ui.input_numeric(
                                "nn_figure_height",
                                "Figure Height",
                                value=6,
                                min=3,
                                max=15
                            ),
                            ui.input_numeric(
                                "nn_figure_dpi",
                                "Figure DPI",
                                value=300,
                                min=72,
                                max=600
                            ),
                            ui.input_numeric(
                                "nn_font_size",
                                "Font Size",
                                value=12,
                                min=8,
                                max=20
                            ),
                            
                            ui.br(),
                            ui.input_action_button(
                                "go_nn_viz",
                                "Generate Visualization",
                                class_="btn-success"
                            ),
                            ui.div(
                                {"style": "padding-top: 20px;"},
                                ui.output_ui("download_button_ui_nn")
                            ),
                        )
                    ),
                    ui.column(
                        9,
                        ui.output_plot(
                            "nn_visualization_plot",
                            width="100%",
                            height="80vh"
                        )
                    )
                ),
            ),
        )
    )

