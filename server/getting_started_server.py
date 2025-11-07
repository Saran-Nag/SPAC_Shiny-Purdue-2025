from shiny import render, reactive, ui

def getting_started_server(input, output, session, shared):
    """Server logic for the Getting Started tab"""
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