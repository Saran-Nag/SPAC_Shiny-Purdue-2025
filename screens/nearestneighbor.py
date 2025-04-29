from shiny import ui

def nearest_neighbor_panel():
    # 10. NEAREST NEIGHBORS PANEL -----------------------------
    return ui.nav_panel("Nearest Neigbors",
        ui.row(
            ui.column(6,
                ui.card(
                    ui.column(12,
                        ui.input_select("nn_anno", "Select an Annotation", choices=[]),
                        ui.output_plot("spac_nearest_neighbor"),
                        ui.output_text_verbatim("nearest_neighbor_profile"),
                        ui.input_action_button("go_nn", "Render Plot", class_="btn-success")
                    )
                ),
            ),
        )
    )
