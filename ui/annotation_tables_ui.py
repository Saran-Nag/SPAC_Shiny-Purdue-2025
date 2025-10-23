"""
UI helper for rendering annotation label tables.

This module contains a single function to build the responsive
annotation tables used by the effect update server. Moving the
markup and styling here keeps server code focused on reactivity
and data handling.
"""

from shiny import ui


def build_annotation_tables(annotation_counts):
    """Build responsive annotation tables markup.

    Parameters
    ----------
    annotation_counts : dict
        Mapping of annotation name -> dict(label -> count).

    Returns
    -------
    shiny.ui.Tag
        A UI element containing the tables and inline styles.
    """
    container = []
    for annotation_name, label_counts_dict in annotation_counts.items():
        # Sort labels by count (descending) and take top 10
        sorted_label_counts = sorted(
            label_counts_dict.items(),
            key=lambda x: x[1],
            reverse=True,
        )[:10]

        # Build table rows
        table_rows = []
        for label, count in sorted_label_counts:
            table_rows.append(
                ui.tags.tr(
                    ui.tags.td(label, {"class": "label-name"}),
                    ui.tags.td(f"{count:,}", {"class": "text-end fw-bold"}),
                )
            )

        # Create a compact table for this annotation
        annotation_table = ui.div(
            {"class": "col-lg-4 col-md-6 col-sm-12 mb-3"},
            ui.div(
                {"class": "annotation-table-card h-100"},
                ui.div(
                    {"class": "table-header"},
                    ui.h6(
                        annotation_name,
                        {"class": "mb-0 text-primary fw-bold"},
                    ),
                ),
                ui.div(
                    {"class": "table-responsive"},
                    ui.tags.table(
                        {"class": "table table-sm table-hover mb-0"},
                        ui.tags.thead(
                            ui.tags.tr(
                                ui.tags.th(
                                    "Label",
                                    {"class": "label-header"},
                                ),
                                ui.tags.th(
                                    "Cells",
                                    {"class": "text-end count-header"},
                                ),
                            ),
                        ),
                        ui.tags.tbody(*table_rows),
                    ),
                ),
            ),
        )
        container.append(annotation_table)

    return ui.div(
        {"class": "row annotation-tables"},
        *container,
        ui.tags.style("""
            .annotation-table-card {
                border: 1px solid #e9ecef;
                border-radius: 0.5rem;
                padding: 0;
                background: white;
                box-shadow: 0 1px 3px rgba(0,0,0,0.1);
                transition: box-shadow 0.2s ease;
            }

            .annotation-table-card:hover {
                box-shadow: 0 4px 8px rgba(0,0,0,0.15);
            }

            .table-header {
                padding: 0.75rem 1rem;
                background: linear-gradient(135deg, #f8f9fa 0%,
                                             #e9ecef 100%);
                border-bottom: 1px solid #dee2e6;
                border-radius: 0.5rem 0.5rem 0 0;
            }

            .annotation-tables .table {
                font-size: 0.85rem;
                margin-bottom: 0;
            }

            .annotation-tables .table th {
                background-color: #f8f9fa;
                border-bottom: 2px solid #dee2e6;
                font-weight: 600;
                font-size: 0.8rem;
                padding: 0.5rem 0.75rem;
            }

            .annotation-tables .table td {
                padding: 0.4rem 0.75rem;
                vertical-align: middle;
            }

            .label-name {
                font-weight: 500;
                color: #495057;
            }

            .annotation-tables .table tbody tr:hover {
                background-color: #f1f3f4;
            }

            .label-header {
                color: #6c757d;
            }

            .count-header {
                color: #6c757d;
            }

            @media (max-width: 768px) {
                .annotation-tables .col-lg-4 {
                    margin-bottom: 1rem;
                }
            }
        """),
    )
