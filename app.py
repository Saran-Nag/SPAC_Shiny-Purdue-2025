from shiny import App, ui, reactive
import os


# Screen imports
from ui import (
    getting_started_ui,
    data_input_ui,
    annotations_ui,
    features_ui,
    boxplot_ui,
    feat_vs_anno_ui,
    anno_vs_anno_ui,
    spatial_ui,
    umap_ui,
    scatterplot_ui,
    nearest_neighbor_ui,
    ripleyL_ui
)   

# Server imports
from server import (
    getting_started_server,
    data_input_server,
    effect_update_server,
    annotations_server,
    features_server,
    boxplot_server,
    feat_vs_anno_server,
    anno_vs_anno_server,
    spatial_server,
    umap_server,
    scatterplot_server,
    nearest_neighbor_server,
    ripleyL_server
)

# Util imports
from utils.data_processing import load_data, read_html_file
from utils.accessibility import accessible_navigation, apply_slider_accessibility_global
from utils.security import apply_security_enhancements


# Read header and footer HTML content
header_html = read_html_file("header.html")
footer_html = read_html_file("footer.html")

file_path = "dev_example.pickle"  # Path to your preloaded .pickle file
preloaded_data = load_data(file_path)  # Initialize as None


app_ui = ui.page_fluid(
    # Apply security enhancements
    apply_security_enhancements(),
    
    # Include header HTML
    ui.HTML(header_html),
    
    # Add navigation accessibility fixes
    accessible_navigation(),
    
    # Add global slider accessibility fixes
    apply_slider_accessibility_global(),
    
    # Main application content
    ui.navset_card_tab(
        getting_started_ui(),
        data_input_ui(),
        annotations_ui(),
        features_ui(),
        boxplot_ui(),
        feat_vs_anno_ui(),
        anno_vs_anno_ui(),
        spatial_ui(),
        umap_ui(),
        scatterplot_ui(),
        nearest_neighbor_ui(),
        ripleyL_ui(),
    ),
    ui.input_action_button(
        "my_fixed_btn",
        "💬",
        #alt = ui.tags.img(
            #src="path/to/your/image.png",
            src="💬",
        #),
        class_="fixed-button btn btn-primary",
    ),
    ui.tags.style("""
        .fixed-button {
            position: fixed;
            bottom: 20px;
            right: 15px;
            z-index: 1000;
            border-radius: 50%;
            width: 60px;
            height: 60px;
            font-size: 15px;
            display: flex;
            align-items: center;
            justify-content: center;
            background: linear-gradient(45deg, #17a2b8, #20c997);
            transition: all 0.3s ease;
            border: none;
            outline: none;
        }
        .fixed-button:hover {
            transform: scale(1.3);
        }
    """),
    # Include footer HTML
    ui.HTML(footer_html)
)


def server(input, output, session):
    # Handle the fixed button click
    @reactive.effect
    @reactive.event(input.my_fixed_btn)
    def show_input_modal():
        # Show a modal dialog when button is clicked
        m = ui.modal(
            ui.h3("Chat with SPAC!"),
            ui.input_text_area(
                "user_input",
                "Type your message here:",
                placeholder="Enter text...",
                rows=8,
                width="100%"
            ),
            easy_close=True,
            footer=ui.div(
                ui.input_action_button("submit_input", "Submit", class_="btn-primary"),
                ui.modal_button("Cancel")
            )
        )
        ui.modal_show(m)

    # Handle the submit button in the modal
    @reactive.effect
    @reactive.event(input.submit_input)
    def handle_submission():
        user_text = input.user_input()
        print(f"Submitted: {user_text}")

        # Close the modal
        ui.modal_remove()
        # Show a confirmation notification
        if user_text is not "":
            ui.notification_show(
                f"Submitted: {user_text}",
                type="message",
                duration=3
            )
        else:
            ui.notification_show(
                f"No Input",
                type="message",
                duration=3
            )
    # Define a reactive variable to track if data is loaded
    data_loaded = reactive.Value(False)

    # Create a reactive variable for the main data
    adata_main = reactive.Value(preloaded_data)  # Initialize with preloaded data

    data_keys = [
        "X_data",
        "obs_data",  # AKA Annotations
        "obsm_data",
        "layers_data",
        "var_data",  # AKA Features
        "uns_data",
        "shape_data",
        "obs_names",
        "obsm_names",
        "layers_names",
        "var_names",
        "uns_names",
        "spatial_distance_columns",
        "df_heatmap",
        "df_relational",
        "df_boxplot",
        "df_histogram2",
        "df_histogram1",
        "df_nn",
        "df_ripley"
    ]

    shared = {
        "preloaded_data": preloaded_data,  # Preloaded data for initial load
        "data_loaded": data_loaded,  # Reactive to track if data is loaded
        "adata_main": adata_main,  # Main anndata object
    }

    # Dynamically create the reactive values for parts of the anndata object
    # and add them to the shared dictionary
    for key in data_keys:
        shared[key] = reactive.Value(None)

    # Individual server components
    getting_started_server(input, output, session, shared)
    
    data_input_server(input, output, session, shared)

    effect_update_server(input, output, session, shared)

    annotations_server(input, output, session, shared)

    features_server(input, output, session, shared)

    boxplot_server(input, output, session, shared)

    feat_vs_anno_server(input, output, session, shared)

    anno_vs_anno_server(input, output, session, shared)

    spatial_server(input, output, session, shared)

    umap_server(input, output, session, shared)

    scatterplot_server(input, output, session, shared)

    nearest_neighbor_server(input, output, session, shared)

    ripleyL_server(input, output, session, shared)


# Create the app with static file serving for www directory
static_path = os.path.join(os.path.dirname(__file__), "www")
app = App(app_ui, server, static_assets=static_path)
