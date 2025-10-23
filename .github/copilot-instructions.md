# Copilot Instructions for SPAC Shiny for Python Development

## Project Overview
This project is a Shiny for Python application for analysis of single cell dataset. Follow these guidelines to maintain code quality, accessibility, and modularity.

## Core Development Principles

### 1. Modular Architecture - UI and Server Separation
- **ALWAYS** separate UI components from server logic
- Place UI modules in the `ui/` directory (e.g., `ui/spatial_ui.py`)
- Place server logic in the `server/` directory (e.g., `server/spatial_server.py`)
- Each feature should have its own dedicated UI and server module
- Keep modules focused on a single responsibility
- Make code easily maintainable by data scientists with clear, readable functions

Example structure:
```python
# ui/feature_ui.py
from shiny import ui

def feature_ui():
    return ui.div(
        # UI components here
    )

# server/feature_server.py  
def feature_server(input, output, session):
    # Server logic here
    pass
```

### 2. Utility Functions and Code Reuse
- **ALWAYS** check the `utils/` folder for existing utilities before creating new functions
- Create reusable functions in appropriate utility modules:
  - `utils/data_processing.py` - Data manipulation and analysis functions
  - `utils/styling.py` - CSS and styling utilities
  - `utils/accessibility.py` - Accessibility helper functions
- When creating new utilities, ensure they are:
  - Well-documented with NumPy-style docstrings
  - Generic enough for reuse across multiple modules
  - Unit testable

### 3. Accessibility Requirements
- **ALWAYS** ensure proper color contrast ratios (WCAG AA: 4.5:1 for normal text, 3:1 for large text)
- Make all interactive elements focusable and keyboard accessible
- Use semantic HTML elements and proper ARIA labels
- Provide alternative text for images and charts
- Use the utilities in `utils/accessibility.py` for consistent accessibility features
- Test with screen readers and keyboard-only navigation

Example accessibility practices:
```python
# Good - accessible button with proper contrast and focus
ui.input_action_button(
    "btn_analyze", 
    "Analyze Data",
    class_="btn-primary accessible-button",
    **{"aria-label": "Start data analysis process"}
)

# Good - accessible plot with alt text
ui.output_plot("plot", alt="Scatter plot showing cell distribution")
```

### 4. Documentation Standards - NumPy Style
- Use NumPy-style docstrings for all functions, classes, and modules
- Include clear parameter descriptions, return values, and examples

```python
def process_spatial_data(data, method="default", threshold=0.5):
    """
    Process spatial data using specified method and threshold.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Input spatial data with x, y coordinates and cell features
    method : str, optional
        Processing method to use, by default "default"
    threshold : float, optional
        Threshold value for filtering, by default 0.5
        
    Returns
    -------
    pandas.DataFrame
        Processed spatial data with additional calculated columns
        
    Raises
    ------
    ValueError
        If data is empty or missing required columns
        
    Examples
    --------
    >>> df = pd.DataFrame({"x": [1, 2], "y": [3, 4], "feature": [0.6, 0.8]})
    >>> result = process_spatial_data(df, threshold=0.7)
    >>> len(result)
    1
    """
```

### 5. Code Style - PEP 8 Compliance
- Follow PEP 8 for all Python code
- Use 4 spaces for indentation
- Line length limit: 79 characters for code, 72 for docstrings
- Use meaningful variable and function names
- Import organization: standard library → third-party → local imports

```python
# Good PEP 8 example
import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
from shiny import App, ui, render

from utils.data_processing import load_spatial_data
from utils.accessibility import create_accessible_plot
```

### 6. Commit Message Format
Use conventional commits with this exact format:
```
<type>(scope): <short description>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, no logic change)
- `refactor`: Code refactoring
- `test`: Adding or modifying tests
- `chore`: Maintenance tasks

**Examples:**
```
feat(spatial): add new nearest neighbor analysis module
fix(ui): correct color contrast issues in plot legends
docs(api): update docstrings for data processing functions
refactor(server): separate plot generation into utility function
style(ui): format code according to PEP 8 standards
test(utils): add unit tests for accessibility functions
```

## File Organization Guidelines

### Directory Structure
```
/
├── ui/                 # All UI modules
├── server/            # All server logic modules  
├── utils/             # Shared utility functions
│   ├── data_processing.py
│   ├── styling.py
│   └── accessibility.py
├── www/               # Static assets (images, CSS, JS)
├── app.py             # Main application entry point
└── requirements.txt   # Dependencies
```

### Naming Conventions
- Use snake_case for files, functions, and variables
- Use descriptive names that indicate functionality
- UI modules: `feature_name_ui.py`
- Server modules: `feature_name_server.py`
- Utility modules: `functionality_name.py`

## Testing and Quality Assurance

### Code Quality Checklist
- [ ] All functions have NumPy-style docstrings
- [ ] Code follows PEP 8 standards
- [ ] UI and server logic are properly separated
- [ ] Accessibility requirements are met
- [ ] Reusable functions are placed in utils/
- [ ] Commit messages follow conventional format
- [ ] No hardcoded values (use configuration)
- [ ] Error handling is implemented
- [ ] Code is readable by data scientists

### Performance Considerations
- Use lazy loading for large datasets
- Implement caching for expensive computations
- Optimize plot rendering for large datasets
- Consider using reactive calculations efficiently

## Shiny for Python Specific Guidelines

### Reactive Programming
- Use `@reactive.calc` for expensive computations
- Use `@reactive.effect` for side effects
- Minimize reactive dependencies
- Use `req()` for input validation

### UI Best Practices
- Use consistent spacing and layout
- Implement responsive design
- Provide loading indicators for long operations
- Use appropriate input types for data

### Server Best Practices
- Validate inputs before processing
- Handle errors gracefully with user-friendly messages
- Use session state appropriately
- Implement proper data flow patterns

## Example Module Template

```python
# ui/example_ui.py
"""
Example UI module for demonstration purposes.

This module provides the user interface components for example functionality.
"""

from shiny import ui
from utils.accessibility import create_accessible_input

def example_ui():
    """
    Create the example feature UI.
    
    Returns
    -------
    shiny.ui.TagList
        UI components for the example feature
    """
    return ui.div(
        ui.h3("Example Feature", class_="accessible-heading"),
        create_accessible_input("example_input", "Enter value:"),
        ui.output_plot("example_plot")
    )

# server/example_server.py
"""
Example server module for demonstration purposes.

This module handles the server-side logic for example functionality.
"""

import pandas as pd
from shiny import reactive, render, req
from utils.data_processing import process_example_data
from utils.accessibility import create_accessible_plot

def example_server(input, output, session):
    """
    Server logic for example feature.
    
    Parameters
    ----------
    input : shiny.session.Inputs
        Shiny input object
    output : shiny.session.Outputs  
        Shiny output object
    session : shiny.session.Session
        Shiny session object
    """
    
    @reactive.calc
    def processed_data():
        """Process input data reactively."""
        req(input.example_input())
        return process_example_data(input.example_input())
    
    @output
    @render.plot
    def example_plot():
        """Render accessible example plot."""
        data = processed_data()
        return create_accessible_plot(data, title="Example Analysis")
```

Remember: The goal is to create maintainable, accessible, and well-documented code that data scientists can easily understand and extend.