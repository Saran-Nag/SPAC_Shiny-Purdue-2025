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
    # Fixed chat button
    ui.input_action_button(
        "my_fixed_btn",
        "💬",
        class_="fixed-button btn btn-primary",
        onclick="toggleChatPanel()"  # Direct JavaScript call
    ),
    # Floating chat panel
    ui.div(
        ui.div(
            ui.h3("Chat with SPAC!", style="margin-top: 0;"),
            ui.tags.button(
                "×",
                onclick="document.getElementById('chat_panel').style.display='none'",
                class_="close-chat-btn",
                type="button"
            ),
            style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;"
        ),
        ui.input_text_area(
            "user_input",
            "Type your message here:",
            placeholder="Enter text...",
            rows=8,
            width="100%"
        ),
        ui.input_action_button("submit_input", "Submit", class_="btn-primary", style="width: 100%; margin-top: 10px;"),
        id="chat_panel",
        class_="chat-panel",
        style="display: none;"  # Hidden by default
    ),
    # JavaScript for toggling
    ui.tags.script("""
        function toggleChatPanel() {
            var panel = document.getElementById('chat_panel');
            if (panel.style.display === 'none' || panel.style.display === '') {
                panel.style.display = 'block';
            } else {
                panel.style.display = 'none';
            }
        }
    """),
    # Styles
    ui.tags.style("""
        .fixed-button {
            position: fixed;
            bottom: 20px;
            right: 15px;
            z-index: 1000;
            border-radius: 50%;
            width: 60px;
            height: 60px;
            min-width: 60px;
            min-height: 60px;
            max-width: 60px;
            max-height: 60px;
            font-size: 24px;
            text-align: center;
            line-height: 60px;
            padding: 0;
            background: linear-gradient(45deg, #17a2b8, #20c997);
            transition: all 0.3s ease;
            border: none;
            outline: none;
            box-sizing: border-box;
        }
        .fixed-button:hover {
            transform: scale(1.2);
            outline: none !important;
            box-shadow: none !important;
        }
        .fixed-button:focus {
            outline: none !important;
            box-shadow: none !important;
        }
        .chat-panel {
            position: fixed;
            bottom: 90px;
            right: 15px;
            width: 350px;
            max-width: 90vw;
            background: white;
            border-radius: 10px;
            padding: 20px;
            z-index: 999;
        }
        .close-chat-btn {
            background: none;
            border: none;
            font-size: 24px;
            cursor: pointer;
            padding: 0;
            width: 30px;
            height: 30px;
            color: #666666;
        }
        .close-chat-btn:hover {
            color: #000000;
            transform: scale(1.2);
        }
    """),
    # Include footer HTML
    ui.HTML(footer_html)
)

def server(input, output, session):
    # Handle the submit button
    @reactive.effect
    @reactive.event(input.submit_input)
    def handle_submission():
        user_text = input.user_input()
        print(f"Submitted: {user_text}")

        # Show notification
        if user_text != "":
            ui.notification_show(
                f"Submitted: {user_text}",
                type="message",
                duration=3
            )
            # Clear the input after submission
            ui.update_text_area("user_input", value="")
        else:
            ui.notification_show(
                "No Input",
                type="warning",
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