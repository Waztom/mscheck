import logging
import os
import sys
import traceback
from logging.handlers import RotatingFileHandler


# Define ANSI color codes for colored console output
class Colors:
    RESET = "\033[0m"
    BOLD = "\033[1m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    BRIGHT_WHITE = "\033[97m"


# Define log format
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

# Define log levels with their names for easier reference
LOG_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


# Create a custom formatter with colors
class ColoredFormatter(logging.Formatter):
    """Custom formatter for colored console output"""

    # Define colors for different log levels
    LEVEL_COLORS = {
        logging.DEBUG: Colors.BRIGHT_BLUE,
        logging.INFO: Colors.BRIGHT_GREEN,
        logging.WARNING: Colors.BRIGHT_YELLOW,
        logging.ERROR: Colors.BRIGHT_RED,
        logging.CRITICAL: Colors.RED + Colors.BOLD,
    }

    def format(self, record):
        # Save original levelname to restore it later
        orig_levelname = record.levelname

        # Add color to the levelname
        levelname_color = self.LEVEL_COLORS.get(record.levelno, Colors.RESET)
        record.levelname = f"{levelname_color}{record.levelname}{Colors.RESET}"

        # Format the message
        result = super().format(record)

        # Restore original levelname
        record.levelname = orig_levelname

        # Add color to specific parts of the message based on content
        if "Analyzing" in result:
            result = result.replace(
                "Analyzing", f"{Colors.BLUE}Analyzing{Colors.RESET}"
            )

        if "Processing sample" in result:
            result = result.replace(
                "Processing sample", f"{Colors.CYAN}Processing sample{Colors.RESET}"
            )

        if "EIC Area:" in result:
            # Highlight high EIC values in green, low in yellow
            parts = result.split("EIC Area: ")
            if len(parts) > 1:
                try:
                    value = float(parts[1].strip())
                    if value > 50.0:
                        result = f"{parts[0]}EIC Area: {Colors.BRIGHT_GREEN}{value}{Colors.RESET}"
                    elif value > 5.0:
                        result = (
                            f"{parts[0]}EIC Area: {Colors.GREEN}{value}{Colors.RESET}"
                        )
                    else:
                        result = (
                            f"{parts[0]}EIC Area: {Colors.YELLOW}{value}{Colors.RESET}"
                        )
                except (ValueError, IndexError):
                    pass

        if "Successfully" in result:
            result = result.replace(
                "Successfully", f"{Colors.BRIGHT_GREEN}Successfully{Colors.RESET}"
            )

        if "Error" in result:
            result = result.replace("Error", f"{Colors.RED}Error{Colors.RESET}")

        if "Warning" in result:
            result = result.replace("Warning", f"{Colors.YELLOW}Warning{Colors.RESET}")

        # Highlight compound types
        for compound_type in ["product", "reactant", "internal-std", "intermediate"]:
            if compound_type in result:
                color = {
                    "product": Colors.MAGENTA,
                    "reactant": Colors.CYAN,
                    "internal-std": Colors.BRIGHT_BLUE,
                    "intermediate": Colors.YELLOW,
                }.get(compound_type, Colors.RESET)
                result = result.replace(
                    compound_type, f"{color}{compound_type}{Colors.RESET}"
                )

        return result


# Add a utility function to log exceptions with traceback
def log_exception(logger, msg="An exception occurred", exc_info=None):
    """
    Log an exception with full traceback information

    Args:
        logger: Logger instance to use
        msg: Error message prefix
        exc_info: Exception info from sys.exc_info() (if None, gets current exception)
    """
    if exc_info is None:
        exc_info = sys.exc_info()

    if exc_info[0] is not None:  # If there's an actual exception
        exc_type, exc_value, exc_traceback = exc_info

        # Log the error message with exception info
        logger.error(f"{msg}: {str(exc_value)}")

        # Format the traceback and log each line
        tb_lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        tb_text = "".join(tb_lines)

        # Log the full traceback at DEBUG level
        logger.debug(f"Traceback:\n{tb_text}")

        # Extract the most relevant traceback info for ERROR level
        tb_summary = traceback.extract_tb(exc_traceback)
        if tb_summary:
            # Get the most recent frame (where the error happened)
            frame = tb_summary[-1]
            file_name = os.path.basename(frame.filename)
            line_number = frame.lineno
            function = frame.name
            code = frame.line

            error_location = f'File "{file_name}", line {line_number}, in {function}'
            logger.error(f"Error occurred at: {error_location}")
            if code:
                logger.error(f"Code: {code}")
    else:
        # No exception context, just log the message
        logger.error(msg)


# Store the root logger configuration
ROOT_LOGGER_CONFIGURED = False


def setup_logger(name="mscheck", level="INFO", log_file=None, use_colors=True):
    """
    Configure the root logger that will be inherited by all loggers
    """
    global ROOT_LOGGER_CONFIGURED

    # Convert string level to logging constant
    numeric_level = getattr(logging, level.upper(), logging.INFO)

    # Configure the root logger only once
    if not ROOT_LOGGER_CONFIGURED:
        # Configure root logger - all other loggers inherit from this
        root_logger = logging.getLogger()
        root_logger.setLevel(numeric_level)

        # Clear any existing handlers to prevent duplicates
        if root_logger.handlers:
            for handler in root_logger.handlers:
                root_logger.removeHandler(handler)

        # Add console handler with color formatting
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(numeric_level)

        if use_colors:
            console_formatter = ColoredFormatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
        else:
            console_formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )

        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)

        # Add file handler if specified
        if log_file:
            # Make sure directory exists
            os.makedirs(os.path.dirname(log_file), exist_ok=True)

            # Set up file handler
            file_handler = RotatingFileHandler(
                log_file, maxBytes=10 * 1024 * 1024, backupCount=5
            )
            file_handler.setLevel(numeric_level)
            file_formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            file_handler.setFormatter(file_formatter)
            root_logger.addHandler(file_handler)

        ROOT_LOGGER_CONFIGURED = True

    # Return the named logger (will inherit config from root)
    return logging.getLogger(name)


def get_logger(name):
    """
    Get a logger with the given name, creating it if necessary.
    Will inherit configuration from the root logger if setup_logger has been called.
    """
    logger = logging.getLogger(name)

    # If root logger hasn't been configured yet, set a basic console handler
    if not ROOT_LOGGER_CONFIGURED and not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

    return logger
