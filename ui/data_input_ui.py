from shiny import ui


def data_input_ui():
    # 1. DATA INPUT PANEL -----------------------------------
    return ui.nav_panel("Data Input",
        # Enhanced CSS for modern, appealing design
        ui.tags.head(ui.tags.style("""
            .shiny-file-input-progress {
                height: 30px !important; 
                line-height: 30px !important;
            }
            .progress-bar {
                height: 30px !important;
                line-height: 30px !important;
                font-size: 16px !important;
            }
            .progress-bar span {
                white-space: nowrap;
                overflow: visible;
            }
            .metric-output {
                font-size: 18px;
                font-weight: 500;
                color: #495057;
            }
            .data-input-hero {
                background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
                border: 1px solid #dee2e6;
                border-radius: 0.75rem;
                padding: 2rem;
                margin-bottom: 1.5rem;
                text-align: center;
            }
            .data-input-title {
                font-size: 2.2rem;
                font-weight: 600;
                color: #495057;
                margin-bottom: 0.5rem;
            }
            .metric-card {
                transition: transform 0.2s ease-in-out, box-shadow 0.2s ease-in-out;
                border: 1px solid #e9ecef;
                border-radius: 0.5rem;
                background: #ffffff;
            }
            .metric-card:hover {
                transform: translateY(-2px);
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }
            .upload-area {
                border: 2px dashed #dee2e6;
                border-radius: 0.75rem;
                padding: 2rem;
                background: #f8f9fa;
                transition: all 0.3s ease;
            }
            .upload-area:hover {
                border-color: #6c757d;
                background: #e9ecef;
            }
            .control-card {
                background: #ffffff;
                border: 1px solid #e9ecef;
                border-radius: 0.75rem;
                transition: box-shadow 0.2s ease-in-out;
            }
            .control-card:hover {
                box-shadow: 0 4px 12px rgba(0,0,0,0.08);
            }
            .stats-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 1rem;
                margin: 1.5rem 0;
            }
            .action-buttons {
                display: flex;
                gap: 0.75rem;
                flex-wrap: wrap;
                margin: 1rem 0;
            }
            .summary-card {
                background: #ffffff;
                border: 1px solid #e9ecef;
                border-radius: 0.75rem;
                margin-top: 1.5rem;
            }
        """)),
        ui.div(
            {"class": "container-fluid p-4"},
            
            # Hero Section
            ui.div(
                {"class": "data-input-hero"},
                ui.h1("SPAC Interactive Dashboard", {"class": "data-input-title"}),
                ui.p("Upload and explore your spatial single-cell data", 
                     {"class": "text-muted mb-0"})
            ),
            
            # Main Content Area
            ui.div(
                {"class": "row"},
                
                # Upload and Data Overview Section
                ui.div(
                    {"class": "col-lg-8 mb-4"},
                    ui.div(
                        {"class": "card control-card"},
                        ui.div(
                            {"class": "card-body"},
                            ui.h5("📁 Data Upload", {"class": "card-title mb-3"}),
                            ui.div(
                                {"class": "upload-area"},
                                ui.input_file(
                                    "input_file", 
                                    "Choose a file to upload (.h5ad or .pickle):", 
                                    multiple=False, 
                                    width="100%"
                                )
                            )
                        )
                    ),
                    
                    # Data Metrics Grid
                    ui.div(
                        {"class": "stats-grid"},
                        ui.div(
                            {"class": "card metric-card"},
                            ui.div(
                                {"class": "card-body text-center"},
                                ui.div({"class": "metric-output"},
                                    ui.output_text("print_rows")
                                )
                            )
                        ),
                        ui.div(
                            {"class": "card metric-card"},
                            ui.div(
                                {"class": "card-body text-center"},
                                ui.div({"class": "metric-output"},
                                    ui.output_text("print_columns")
                                )
                            )
                        ),
                        ui.div(
                            {"class": "card metric-card"},
                            ui.div(
                                {"class": "card-body text-center"},
                                ui.div({"class": "metric-output"},
                                    ui.output_text("print_obs_names")
                                )
                            )
                        ),
                        ui.div(
                            {"class": "card metric-card"},
                            ui.div(
                                {"class": "card-body text-center"},
                                ui.div({"class": "metric-output"},
                                    ui.output_text("print_obsm_names")
                                )
                            )
                        ),
                        ui.div(
                            {"class": "card metric-card"},
                            ui.div(
                                {"class": "card-body text-center"},
                                ui.div({"class": "metric-output"},
                                    ui.output_text("print_layers_names")
                                )
                            )
                        ),
                        ui.div(
                            {"class": "card metric-card"},
                            ui.div(
                                {"class": "card-body text-center"},
                                ui.div({"class": "metric-output"},
                                    ui.output_text("print_uns_names")
                                )
                            )
                        )
                    )
                ),
                
                # Data Controls Section
                ui.div(
                    {"class": "col-lg-4 mb-4"},
                    ui.div(
                        {"class": "card control-card"},
                        ui.div(
                            {"class": "card-header bg-light"},
                            ui.h5("🔧 Data Controls", {"class": "mb-0"})
                        ),
                        ui.div(
                            {"class": "card-body"},
                            ui.input_checkbox(
                                "subset_select_check", 
                                "Enable Data Subsetting", 
                                False
                            ),
                            ui.div({"class": "mb-3"}, id="main-subset_anno_dropdown"),
                            ui.div({"class": "mb-3"}, id="main-subset_label_dropdown"),
                            
                            ui.div(
                                {"class": "action-buttons"},
                                ui.input_action_button(
                                    "go_subset", 
                                    "Apply Subset", 
                                    class_="btn btn-success flex-fill"
                                ),
                                ui.input_action_button(
                                    "restore_data", 
                                    "Restore Data", 
                                    class_="btn btn-outline-warning flex-fill"
                                )
                            ),
                            
                            ui.div(
                                {"class": "mt-3 p-3 bg-light rounded"},
                                ui.h6("Subset History", {"class": "mb-2"}),
                                ui.div({"class": "metric-output small"},
                                    ui.output_text("print_subset_history")
                                )
                            )
                        )
                    )
                )
            ),
            
            # Annotation Summary Section
            ui.div(
                {"class": "row"},
                ui.div(
                    {"class": "col-12"},
                    ui.div(
                        {"class": "card summary-card"},
                        ui.div(
                            {"class": "card-header bg-light"},
                            ui.h5("📊 Annotation Summary with Top 10 Labels", {"class": "mb-0"})
                        ),
                        ui.div(
                            {"class": "card-body"},
                            ui.output_ui("annotation_labels_display")
                        )
                    )
                )
            )
        )
    )
