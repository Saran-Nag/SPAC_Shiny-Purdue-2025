from shiny import ui
from utils.data_processing import read_markdown_file

def getting_started_ui():
    """UI for the Getting Started tab - landing page with tutorial content"""
    
    return ui.nav_panel(
        "Getting Started",
        ui.tags.style("""
            .hero-section {
                background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
                color: #495057;
                padding: 2.5rem 0;
                margin-bottom: 2rem;
                border-radius: 0.5rem;
                border: 1px solid #dee2e6;
            }
            .feature-card {
                transition: transform 0.2s ease-in-out, box-shadow 0.2s ease-in-out;
                border: 1px solid #e9ecef;
                border-radius: 0.75rem;
                overflow: hidden;
            }
            .feature-card:hover {
                transform: translateY(-2px);
                box-shadow: 0 5px 15px rgba(0,0,0,0.08);
            }
            .feature-icon {
                font-size: 2.5rem;
                margin-bottom: 1rem;
            }
            .workflow-step {
                background: #ffffff;
                border-left: 4px solid #6c757d;
                border: 1px solid #e9ecef;
                margin-bottom: 1rem;
                padding: 1rem;
                border-radius: 0.5rem;
            }
            .stats-card {
                background: linear-gradient(45deg, #6c757d, #868e96);
                color: white;
                text-align: center;
                padding: 1.5rem;
                border-radius: 0.75rem;
            }
            .stats-card-alt {
                background: linear-gradient(45deg, #17a2b8, #20c997);
                color: white;
                text-align: center;
                padding: 1.5rem;
                border-radius: 0.75rem;
            }
            .quick-start-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                gap: 1.5rem;
                margin: 2rem 0;
            }
            .gentle-title {
                font-weight: 600 !important;
                font-size: 2.5rem !important;
                color: #495057 !important;
            }
            .gentle-subtitle {
                font-weight: 400 !important;
                color: #6c757d !important;
            }
        """),
        
        ui.div(
            {"class": "container-fluid p-4"},
            
            # Hero Section
            ui.div(
                {"class": "hero-section text-center"},
                ui.div(
                    {"class": "container"},
                    ui.h1("Welcome to SPAC", {"class": "gentle-title mb-3"}),
                    ui.h2("Spatial Analysis of Cellular Data", {"class": "gentle-subtitle mb-4"}),
                    ui.p("An interactive dashboard for spatial single-cell analysis with 10x performance improvements", 
                         {"class": "lead mb-4"}),
                    ui.div(
                        {"class": "row justify-content-center"},
                        ui.div(
                            {"class": "col-md-3 col-6 mb-3"},
                            ui.div(
                                {"class": "stats-card-alt"},
                                ui.h3("10x", {"class": "fw-bold mb-0"}),
                                ui.p("Faster Analysis", {"class": "mb-0 small"})
                            )
                        ),
                        ui.div(
                            {"class": "col-md-3 col-6 mb-3"},
                            ui.div(
                                {"class": "stats-card"},
                                ui.h3("1M+", {"class": "fw-bold mb-0"}),
                                ui.p("Cells Supported", {"class": "mb-0 small"})
                            )
                        ),
                        ui.div(
                            {"class": "col-md-3 col-6 mb-3"},
                            ui.div(
                                {"class": "stats-card-alt"},
                                ui.h3("8+", {"class": "fw-bold mb-0"}),
                                ui.p("Analysis Tools", {"class": "mb-0 small"})
                            )
                        )
                    )
                )
            ),
            
            # Quick Start Section
            ui.div(
                {"class": "row mb-5"},
                ui.div(
                    {"class": "col-12"},
                    ui.h2("🚀 Quick Start Guide", {"class": "text-center mb-4"}),
                    ui.div(
                        {"class": "quick-start-grid"},
                        
                        # Step 1
                        ui.div(
                            {"class": "workflow-step"},
                            ui.h5("📁 1. Load Your Data", {"class": "fw-bold text-primary mb-2"}),
                            ui.tags.ul(
                                ui.tags.li("Navigate to the Data Input tab"),
                                ui.tags.li("Upload .h5ad or .pickle files"),
                                ui.tags.li("Sample data pre-loaded for demo"),
                                class_="mb-0"
                            )
                        ),
                        
                        # Step 2
                        ui.div(
                            {"class": "workflow-step"},
                            ui.h5("🔍 2. Explore Annotations", {"class": "fw-bold text-primary mb-2"}),
                            ui.tags.ul(
                                ui.tags.li("Check cell type annotations"),
                                ui.tags.li("Review metadata and conditions"),
                                ui.tags.li("Filter and subset your data"),
                                class_="mb-0"
                            )
                        ),
                        
                        # Step 3
                        ui.div(
                            {"class": "workflow-step"},
                            ui.h5("📊 3. Analyze Features", {"class": "fw-bold text-primary mb-2"}),
                            ui.tags.ul(
                                ui.tags.li("Explore gene expression patterns"),
                                ui.tags.li("Create statistical visualizations"),
                                ui.tags.li("Compare features across conditions"),
                                class_="mb-0"
                            )
                        ),
                        
                        # Step 4
                        ui.div(
                            {"class": "workflow-step"},
                            ui.h5("🗺️ 4. Spatial Analysis", {"class": "fw-bold text-primary mb-2"}),
                            ui.tags.ul(
                                ui.tags.li("View spatial distribution plots"),
                                ui.tags.li("Perform neighborhood analysis"),
                                ui.tags.li("Calculate spatial statistics"),
                                class_="mb-0"
                            )
                        )
                    )
                )
            ),
            
            # Features Section
            ui.div(
                {"class": "row mb-5"},
                ui.div(
                    {"class": "col-12"},
                    ui.h2("✨ Key Features", {"class": "text-center mb-4"}),
                    ui.div(
                        {"class": "row"},
                        
                        # Spatial Analysis
                        ui.div(
                            {"class": "col-lg-3 col-md-6 mb-4"},
                            ui.div(
                                {"class": "card feature-card h-100 shadow-sm"},
                                ui.div(
                                    {"class": "card-body text-center"},
                                    ui.div("🔬", {"class": "feature-icon"}),
                                    ui.h5("Spatial Analysis", {"class": "card-title"}),
                                    ui.p("Interactive spatial plots, neighborhood analysis, and Ripley's L function",
                                         {"class": "card-text"})
                                )
                            )
                        ),
                        
                        # Performance
                        ui.div(
                            {"class": "col-lg-3 col-md-6 mb-4"},
                            ui.div(
                                {"class": "card feature-card h-100 shadow-sm"},
                                ui.div(
                                    {"class": "card-body text-center"},
                                    ui.div("⚡", {"class": "feature-icon"}),
                                    ui.h5("High Performance", {"class": "card-title"}),
                                    ui.p("10x faster boxplots and histograms, optimized for large datasets",
                                         {"class": "card-text"})
                                )
                            )
                        ),
                        
                        # Visualizations
                        ui.div(
                            {"class": "col-lg-3 col-md-6 mb-4"},
                            ui.div(
                                {"class": "card feature-card h-100 shadow-sm"},
                                ui.div(
                                    {"class": "card-body text-center"},
                                    ui.div("📈", {"class": "feature-icon"}),
                                    ui.h5("Rich Visualizations", {"class": "card-title"}),
                                    ui.p("UMAP, scatter plots, heatmaps with interactive features",
                                         {"class": "card-text"})
                                )
                            )
                        ),
                        
                        # Data Formats
                        ui.div(
                            {"class": "col-lg-3 col-md-6 mb-4"},
                            ui.div(
                                {"class": "card feature-card h-100 shadow-sm"},
                                ui.div(
                                    {"class": "card-body text-center"},
                                    ui.div("💾", {"class": "feature-icon"}),
                                    ui.h5("Flexible Input", {"class": "card-title"}),
                                    ui.p("Supports AnnData (.h5ad) and pickle formats with spatial coordinates",
                                         {"class": "card-text"})
                                )
                            )
                        )
                    )
                )
            ),
            
            # SPAC Terminology Section
            ui.div(
                {"class": "row mb-5"},
                ui.div(
                    {"class": "col-12"},
                    ui.h2("📚 SPAC Terminology", {"class": "text-center mb-4"}),
                    ui.div(
                        {"class": "card shadow-sm"},
                        ui.div(
                            {"class": "card-body"},
                            ui.p("SPAC uses general terminology to simplify technical terms from the AnnData object for less technical users. Here is a quick guide:", 
                                 {"class": "lead mb-4"}),
                            ui.div(
                                {"class": "row"},
                                ui.div(
                                    {"class": "col-md-6"},
                                    ui.tags.ul(
                                        ui.tags.li([
                                            ui.tags.strong("Cells: "), 
                                            "Rows in the X matrix of AnnData."
                                        ], {"class": "mb-2"}),
                                        ui.tags.li([
                                            ui.tags.strong("Features: "), 
                                            "Columns in the X matrix of AnnData, representing gene expression or antibody intensity."
                                        ], {"class": "mb-2"}),
                                        ui.tags.li([
                                            ui.tags.strong("Tables: "), 
                                            "Originally called layers in AnnData, these represent transformed features."
                                        ], {"class": "mb-2"}),
                                        {"class": "list-unstyled"}
                                    )
                                ),
                                ui.div(
                                    {"class": "col-md-6"},
                                    ui.tags.ul(
                                        ui.tags.li([
                                            ui.tags.strong("Associated Tables: "), 
                                            "Corresponds to .obsm in AnnData and can store spatial coordinates, UMAP embeddings, etc."
                                        ], {"class": "mb-2"}),
                                        ui.tags.li([
                                            ui.tags.strong("Annotations: "), 
                                            "Corresponds to .obs in AnnData and can store cell phenotypes, experiment names, slide IDs, etc."
                                        ], {"class": "mb-2"}),
                                        {"class": "list-unstyled"}
                                    )
                                )
                            ),
                            ui.div(
                                {"class": "text-center mt-3 pt-3 border-top"},
                                ui.p([
                                    "For more in-depth explanations, visit our ",
                                    ui.a("GitHub page", 
                                         href="https://github.com/FNLCR-DMAP/spac_datamine/blob/main/CONTRIBUTING.md", 
                                         target="_blank",
                                         class_="text-decoration-none"),
                                    "."
                                ], {"class": "mb-0"})
                            )
                        )
                    )
                )
            ),
            
            # Tips Section
            ui.div(
                {"class": "row"},
                ui.div(
                    {"class": "col-md-6 mb-4"},
                    ui.div(
                        {"class": "card border-info"},
                        ui.div(
                            {"class": "card-header bg-light text-dark border-info"},
                            ui.h5("💡 Pro Tips", {"class": "mb-0"})
                        ),
                        ui.div(
                            {"class": "card-body"},
                            ui.tags.ul(
                                ui.tags.li("Start with the Data Input tab to load your dataset"),
                                ui.tags.li("Use filtering options for large datasets (>1M cells)"),
                                ui.tags.li("Hover over controls for helpful tooltips"),
                                ui.tags.li("Export plots for publication-ready figures"),
                                class_="mb-0"
                            )
                        )
                    )
                ),
                ui.div(
                    {"class": "col-md-6 mb-4"},
                    ui.div(
                        {"class": "card border-secondary"},
                        ui.div(
                            {"class": "card-header bg-light text-dark border-secondary"},
                            ui.h5("📋 Data Requirements", {"class": "mb-0"})
                        ),
                        ui.div(
                            {"class": "card-body"},
                            ui.tags.ul(
                                ui.tags.li("AnnData objects with spatial coordinates in obsm['spatial']"),
                                ui.tags.li("Cell annotations in adata.obs"),
                                ui.tags.li("Gene expression data in adata.X"),
                                ui.tags.li("Supports .h5ad and .pickle file formats"),
                                class_="mb-0"
                            )
                        )
                    )
                )
            ),
            
            # Footer
            ui.div(
                {"class": "text-center mt-5 pt-4 border-top"},
                ui.p("Built with ❤️ by the Frederick National Laboratory for Cancer Research and Purdue University",
                     {"class": "text-muted"})
            )
        )
    )