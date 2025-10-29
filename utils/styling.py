"""
Styling utilities for SPAC Shiny application.

This module provides centralized CSS styles and styling functions
to maintain consistency and improve code organization.
"""

from shiny import ui


def get_data_input_styles():
    """
    Get CSS styles specific to the Data Input page.
    
    Returns:
        ui.tags.style: CSS styles for data input components
    """
    return ui.tags.style("""
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
            border-radius: 1rem;
            padding: 2.5rem 2rem;
            margin-bottom: 1.5rem;
            text-align: center;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
        }
        .data-input-title {
            font-size: 1.8rem;
            font-weight: 500;
            color: #343a40;
            margin-bottom: 0.75rem;
            letter-spacing: -0.02em;
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 0.5rem;
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
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
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
        .metric-list {
            background: #f8f9fa;
            border-radius: 0.375rem;
            padding: 0.75rem;
            font-family: 'Segoe UI', 'Arial', sans-serif;
            font-size: 0.8rem;
            line-height: 1.3;
            color: #212529;
            max-height: 120px;
            overflow-y: auto;
        }
        .metric-item {
            display: flex;
            align-items: flex-start;
            margin-bottom: 0.25rem;
            word-break: break-word;
            hyphens: auto;
        }
        .metric-bullet {
            color: #212529;
            margin-right: 0.5rem;
            flex-shrink: 0;
            line-height: 1.3;
        }
        .metric-text {
            flex: 1;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .metric-title {
            font-size: 0.9rem !important;
            font-weight: 600 !important;
            margin-bottom: 0.5rem !important;
        }
        
        /* Specific accessibility improvements for data input cards */
        .metric-card h6 {
            color: #212529 !important;
        }
        
        .metric-card .text-success,
        .metric-card .fw-semibold.text-success {
            color: #146c43 !important;
        }
        
        .metric-card .text-info,
        .metric-card .fw-semibold.text-info {
            color: #055160 !important;
        }
        
        .metric-card .text-warning,
        .metric-card .fw-semibold.text-warning {
            color: #664d03 !important;
        }
        
        .metric-card .text-secondary,
        .metric-card .fw-semibold.text-secondary {
            color: #495057 !important;
        }
    """)


def get_global_styles():
    """
    Get global CSS styles that apply across the entire application.
    
    Returns:
        ui.tags.style: Global CSS styles
    """
    return ui.tags.style("""
        /* Global color scheme and typography */
        :root {
            --primary-color: #0d6efd;
            --secondary-color: #495057;
            --success-color: #146c43;
            --info-color: #055160;
            --warning-color: #664d03;
            --danger-color: #b02a37;
            --light-color: #f8f9fa;
            --dark-color: #212529;
        }
        
        /* Global accessibility improvements */
        .nav-link:focus,
        .btn:focus,
        input:focus,
        select:focus,
        textarea:focus {
            outline: 2px solid var(--primary-color) !important;
            outline-offset: 2px !important;
        }
        
        /* Ensure sufficient color contrast */
        .text-muted {
            color: #495057 !important;
        }
        
        /* Accessible text colors for cards - using dark colors for AA compliance */
        .text-success {
            color: #146c43 !important;
        }
        
        .text-info {
            color: #055160 !important;
        }
        
        .text-warning {
            color: #664d03 !important;
        }
        
        .text-secondary {
            color: #495057 !important;
        }
        
        /* Ensure all card text has good contrast */
        .card-body h6 {
            color: #212529 !important;
        }
        
        /* Button color improvements for accessibility */
        .btn-outline-warning {
            color: #664d03 !important;
            border-color: #664d03 !important;
            background-color: transparent !important;
        }
        
        .btn-outline-warning:hover,
        .btn-outline-warning:focus {
            background-color: #664d03 !important;
            border-color: #664d03 !important;
            color: #ffffff !important;
        }
        
        .btn-outline-warning:active {
            background-color: #664d03 !important;
            border-color: #664d03 !important;
            color: #ffffff !important;
        }
        
        /* Responsive design improvements */
        @media (max-width: 768px) {
            .container-fluid {
                padding: 1rem !important;
            }
            
            .stats-grid {
                grid-template-columns: 1fr !important;
            }
            
            .action-buttons {
                flex-direction: column;
            }
        }
    """)


def get_component_styles():
    """
    Get reusable component styles that can be used across different pages.
    
    Returns:
        ui.tags.style: Component CSS styles
    """
    return ui.tags.style("""
        /* Card components */
        .card-hover {
            transition: transform 0.2s ease-in-out, box-shadow 0.2s ease-in-out;
        }
        
        .card-hover:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }
        
        /* Button groups */
        .btn-group-custom {
            display: flex;
            gap: 0.5rem;
            flex-wrap: wrap;
        }
        
        .btn-group-custom .btn {
            flex: 1;
            min-width: 120px;
        }
        
        /* Loading states */
        .loading-overlay {
            position: relative;
        }
        
        .loading-overlay::after {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(255, 255, 255, 0.8);
            display: flex;
            align-items: center;
            justify-content: center;
        }
        
        /* Form improvements */
        .form-group label {
            font-weight: 500;
            margin-bottom: 0.5rem;
        }
        
        .form-control, .form-select {
            border-radius: 0.375rem;
            border: 1px solid #ced4da;
            transition: border-color 0.2s ease-in-out, box-shadow 0.2s ease-in-out;
        }
        
        .form-control:focus, .form-select:focus {
            border-color: var(--primary-color);
            box-shadow: 0 0 0 0.2rem rgba(13, 110, 253, 0.25);
        }
        
        /* Status indicators */
        .status-indicator {
            display: inline-block;
            width: 8px;
            height: 8px;
            border-radius: 50%;
            margin-right: 0.5rem;
        }
        
        .status-success { background-color: var(--success-color); }
        .status-warning { background-color: var(--warning-color); }
        .status-danger { background-color: var(--danger-color); }
        .status-info { background-color: var(--info-color); }
    """)


def create_styled_head(*additional_styles):
    """
    Create a comprehensive head section with all styling.
    
    Args:
        *additional_styles: Additional ui.tags.style objects to include
        
    Returns:
        ui.tags.head: Complete head section with all styles
    """
    styles = [
        get_global_styles(),
        get_component_styles()
    ]
    
    # Add any additional styles passed in
    styles.extend(additional_styles)
    
    return ui.tags.head(*styles)


def create_data_input_head():
    """
    Create head section specifically for data input page.
    
    Returns:
        ui.tags.head: Head section with data input styles
    """
    return create_styled_head(get_data_input_styles())


# Style presets for common patterns
def hero_section_class():
    """Return CSS class name for hero sections."""
    return "data-input-hero"


def metric_card_class():
    """Return CSS class name for metric cards."""
    return "card metric-card"


def control_card_class():
    """Return CSS class name for control cards."""
    return "card control-card"


def action_buttons_class():
    """Return CSS class name for action button containers."""
    return "action-buttons"


def stats_grid_class():
    """Return CSS class name for statistics grids."""
    return "stats-grid"


def metric_list_class():
    """Return CSS class name for metric lists."""
    return "metric-list"