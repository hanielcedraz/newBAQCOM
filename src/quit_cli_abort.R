quit_cli_abort <- function(message) {
    tryCatch(
        {
            rlang::with_options(
                list(rlang_backtrace_on_error = "none"),
                stop(message)
            )
        },
        error = function(e) {
            # Construct the custom error message with red "X"
            red_x <- "\033[31m✖\033[0m"
            modified_message <- paste("Error:\n", red_x, message)
            # Print the custom error message
            cat(modified_message, "\n")
            # Stop quietly
            invokeRestart("abort")
        }
    )
}
