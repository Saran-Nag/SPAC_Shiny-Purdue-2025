"""
SPAC Shiny App - Data Processing Module
This module contains functions for loading and processing data for the SPAC Shiny app.
"""

import anndata as ad
import pandas as pd


def read_html_file(filepath):
    """Read HTML file content"""
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            return file.read()
    except FileNotFoundError:
        return ""


def read_markdown_file(filepath):
    """Read Markdown file content and convert to HTML"""
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            content = file.read()
            
            def process_inline_formatting(text):
                """Process bold, italic, and code formatting"""
                # Handle bold text **text**
                import re
                text = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', text)
                # Handle italic text *text*
                text = re.sub(r'\*(.*?)\*', r'<em>\1</em>', text)
                # Handle inline code `code`
                text = re.sub(r'`(.*?)`', r'<code>\1</code>', text)
                return text
            
            lines = content.split('\n')
            html_lines = []
            in_ordered_list = False
            in_unordered_list = False
            
            for line in lines:
                original_line = line
                line = line.strip()
                
                # Handle headers
                if line.startswith('#### '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[5:])
                    html_lines.append(f'<h4>{header_text}</h4>')
                elif line.startswith('### '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[4:])
                    html_lines.append(f'<h3>{header_text}</h3>')
                elif line.startswith('## '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[3:])
                    html_lines.append(f'<h2>{header_text}</h2>')
                elif line.startswith('# '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[2:])
                    html_lines.append(f'<h1>{header_text}</h1>')
                
                # Handle numbered lists (1. 2. 3. etc.)
                elif line and line[0].isdigit() and '. ' in line:
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    if not in_ordered_list:
                        html_lines.append('<ol>')
                        in_ordered_list = True
                    list_text = line.split('. ', 1)[1]
                    formatted_text = process_inline_formatting(list_text)
                    html_lines.append(f'<li>{formatted_text}</li>')
                
                # Handle bullet lists (- or *)
                elif line.startswith('- ') or line.startswith('* '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if not in_unordered_list:
                        html_lines.append('<ul>')
                        in_unordered_list = True
                    list_text = line[2:]
                    formatted_text = process_inline_formatting(list_text)
                    html_lines.append(f'<li>{formatted_text}</li>')
                
                # Handle empty lines
                elif not line:
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    html_lines.append('<br>')
                
                # Handle regular paragraphs
                else:
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    formatted_text = process_inline_formatting(line)
                    html_lines.append(f'<p>{formatted_text}</p>')
            
            # Close any remaining lists
            if in_ordered_list:
                html_lines.append('</ol>')
            if in_unordered_list:
                html_lines.append('</ul>')
            
            return '\n'.join(html_lines)
            
    except FileNotFoundError:
        return "<div class='alert alert-warning'>Content not found.</div>"
    except Exception as e:
        return f"<div class='alert alert-danger'>Error: {str(e)}</div>"

def load_data(file_path):
    """
    Load data from a specified file path. Supports .h5ad and .pickle formats.
    
    Parameters:
        file_path (str): The path to the data file.

    Returns:
        adata: Loaded AnnData object or None if loading fails.
    """
    try:
        with open(file_path, 'rb') as file:
            if file_path.endswith('.h5ad'):
                try:
                    adata = ad.read_h5ad(file)
                    return adata
                except Exception as e:
                    print(f"Error loading .h5ad file: {e}")
                    return None
            elif file_path.endswith('.pickle'):
                try:
                    adata = pd.read_pickle(file)
                    return adata
                except Exception as e:
                    print(f"Error loading .pickle file: {e}")
                    return None
            else:
                print("Unsupported file format. Please provide a .h5ad or .pickle file.")
                return None
    except FileNotFoundError:
        print(
            f"Preloaded data file {file_path} not found.",
            "Proceeding without preloaded data."
        )


def get_annotation_label_counts(adata: ad.AnnData):
    """
    Return a dictionary of every annotation (column in adata.obs),
    where the value is a dict of {label: cell_count}.

    Parameters:
        adata (AnnData): AnnData object containing the data.

    Returns:
        dict: A dictionary where keys are annotation labels and values 
        are dictionaries of {label: cell_count} for each annotation.

    Example structure:
      {
        "cell_type": {"T-cell": 100, "B-cell": 80, ...},
        "condition": {"disease": 120, "healthy": 60, ...},
        ...
      }
    """
    if adata is None or not hasattr(adata, "obs") or adata.obs.empty:
        return {}

    annotation_counts = {}
    for col in adata.obs.columns:
        # value_counts returns a Series of {label: count}
        vc = adata.obs[col].value_counts(dropna=False)
        annotation_counts[col] = vc.to_dict()

    return annotation_counts
