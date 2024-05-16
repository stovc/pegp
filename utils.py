from colorama import Fore, Back, Style


def highlight_text(text):
    """Get string. Return highlighted string using colorama."""
    highlighted_text = Fore.BLACK + Back.WHITE + text + Style.RESET_ALL
    return highlighted_text

