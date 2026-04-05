#Code from https://stackoverflow.com/questions/9321741/printing-to-screen-and-writing-to-a-file-at-the-same-time
#Adapted from there
import programme_settings
import os

programme_settings.read()

from typing import Callable

def print_logger(
        old_print: Callable,
        file_name: str,
) -> Callable:
    """Returns a function which calls `old_print` twice, specifying a `file=` on the second call.

    Arguments:
        old_print: The `print` function to call twice.
        file_name: The name to give the log file.
    """

    def log_print(*args, **kwargs):
        old_print(*args, **kwargs)
        with open(file_name, "a") as log_file:
            old_print(*args, file=log_file, **kwargs)

    return log_print

#create a global log object
fle = programme_settings.settings["DEFAULT"]["print_logger"]

if os.path.exists(fle):
                os.remove(fle)
logger = print_logger(print, fle)

