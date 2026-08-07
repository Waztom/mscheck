"""Command-line entry point for running a bulk MSCheck analysis.

Usage:
    python -m mscheck path/to/config.yaml
    python -m mscheck path/to/config.yaml --log-dir logs --batch-size 20
"""
import argparse
import os
from datetime import datetime

from mscheck.bulkanalyse import BulkAnalyser
from mscheck.logging_config import setup_logger, get_logger


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="mscheck",
        description="Run a bulk MSCheck analysis from a YAML configuration file.",
    )
    parser.add_argument("config", help="Path to the YAML configuration file")
    parser.add_argument(
        "--log-dir",
        default="logs",
        help="Directory to write the timestamped log file to (default: %(default)s)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=20,
        help="Number of samples to process per batch (default: %(default)s)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    os.makedirs(args.log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(args.log_dir, f"mscheck_{timestamp}.log")

    setup_logger(name="MSCheck", level="INFO", log_file=log_file, use_colors=True)
    logger = get_logger("MSCheck")
    logger.info(f"Logging to file: {log_file}")

    analyzer = BulkAnalyser(args.config)
    report_paths = analyzer.run_complete_workflow(batch_size=args.batch_size)

    logger.info(f"Analysis complete! Generated {len(report_paths)} reports")
    return report_paths


if __name__ == "__main__":
    main()
