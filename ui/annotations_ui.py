from shiny import ui
from utils.accessibility import accessible_slider


def annotations_ui():
    # 2. ANNOTATIONS PANEL (Histogram of annotations) --------
    return ui.nav_panel(
        "Annotations",
        ui.card(
            {"style": "width:125%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "h2_anno",
                            "Select an Annotation",
                            choices=[]
                        ),
                        ui.input_checkbox(
                            "h2_group_by_check",
                            "Group By",
                            value=False
                        ),
                        ui.div(id="main-h2_dropdown"),
                        ui.div(id="main-h2_check"),
                        ui.div(id="main-h2_together_drop"),
                        accessible_slider(
                            "anno_slider",
                            "Rotate X-axis Labels (degrees)",
                            min_val=0,
                            max_val=90,
                            value=0,
                            step=1
                        ),
                        ui.div(style="height: 20px;"),
                        ui.input_slider(
                            "annotations_font_size",
                            "Axis Label Font Size",
                            min=3,
                            max=24,
                            value=10
                        ),
                        ui.input_action_button(
                            "go_h2",
                            "Render Plot",
                            class_="btn-success"
                        ),
                        ui.div(
                            {"style": "padding-top: 20px;"},
                            ui.output_ui("download_histogram_button_ui")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px;"},
                            ui.output_plot(
                                "spac_Histogram_2",
                                width="100%",
                                height="80vh"
                            )
                        )
                    )
                )
            )
        )
    )
