#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-05 16:22:13
# @DESCRIPTION: Convert BioSample RunInfo files to CSV format
# @VERSION: v0.0.1

import argparse
import csv
import os
import re
import sys
from pathlib import Path


def parse_runinfo_file(file_path):
    """
    Parse a BioSample RunInfo file into a list of dictionaries.
    Each dictionary represents a sample's information.
    """
    samples = []
    current_sample = None

    with open(file_path, "r") as f:
        content = f.read()

    # Split by sample entries (they start with a number followed by ":")
    sample_blocks = re.split(r"^\d+:\s", content, flags=re.MULTILINE)[1:]

    for block in sample_blocks:
        sample_data = {}

        # Extract sample name
        name_match = re.match(r"Sample\s+(.*?)\n", block)
        if name_match:
            sample_data["sample_name"] = name_match.group(1).strip()

        # Extract identifiers
        identifiers_match = re.search(r"Identifiers:\s+(.*?)\n", block)
        if identifiers_match:
            identifiers_line = identifiers_match.group(1).strip()

            # Extract BioSample
            biosample_match = re.search(r"BioSample:\s+(\w+)", identifiers_line)
            if biosample_match:
                sample_data["biosample_id"] = biosample_match.group(1)

            # Extract SRA
            sra_match = re.search(r"SRA:\s+(\w+)", identifiers_line)
            if sra_match:
                sample_data["sra_id"] = sra_match.group(1)

            # Extract GEO
            geo_match = re.search(r"GEO:\s+(\w+)", identifiers_line)
            if geo_match:
                sample_data["geo_id"] = geo_match.group(1)

        # Extract organism
        organism_match = re.search(r"Organism:\s+(.*?)\n", block)
        if organism_match:
            sample_data["organism"] = organism_match.group(1).strip()

        # Extract attributes
        attributes_section = re.search(
            r"Attributes:(.*?)(?:Accession:|$)", block, re.DOTALL
        )
        if attributes_section:
            attributes_text = attributes_section.group(1).strip()
            attributes_lines = [
                line.strip()
                for line in attributes_text.split("\n")
                if line.strip()
            ]

            for attr_line in attributes_lines:
                # Match attribute name and value
                attr_match = re.match(r'/([^=]+)="([^"]+)"', attr_line)
                if attr_match:
                    attr_name = attr_match.group(1).strip()
                    attr_value = attr_match.group(2).strip()
                    # Convert attribute names to valid column names
                    attr_name = attr_name.replace(" ", "_").lower()
                    sample_data[attr_name] = attr_value

        # Extract accession and ID
        accession_match = re.search(r"Accession:\s+(\w+)\s+ID:\s+(\d+)", block)
        if accession_match:
            sample_data["accession"] = accession_match.group(1)
            sample_data["id"] = accession_match.group(2)

        samples.append(sample_data)

    return samples


def write_to_csv(samples, csv_file):
    """
    Write the samples data to a CSV file.
    """
    if not samples:
        return

    # Get all unique keys across all samples
    fieldnames = set()
    for sample in samples:
        fieldnames.update(sample.keys())

    fieldnames = sorted(list(fieldnames))

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(samples)


def validate_file_path(file_path):
    """
    Validate that the input file exists.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    return file_path


def ensure_output_dir(output_file):
    """
    Ensure the output directory exists.
    """
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)


def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Convert BioSample RunInfo file to CSV format."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input BioSample RunInfo file",
        type=validate_file_path,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output CSV file. If not specified, will use the input filename with .csv extension",
    )

    return parser.parse_args()


def main():
    """
    Main function to process BioSample RunInfo file and output CSV.
    """
    args = parse_args()

    # Set input and output files
    runinfo_file = args.input

    if args.output:
        csv_file = args.output
    else:
        csv_file = str(Path(runinfo_file).with_suffix(".csv"))

    # Ensure output directory exists
    ensure_output_dir(csv_file)

    # Parse RunInfo to structured data
    try:
        samples = parse_runinfo_file(runinfo_file)

        # Write to CSV file
        write_to_csv(samples, csv_file)

        print(f"Successfully converted {runinfo_file} to {csv_file}")
        print(f"Processed {len(samples)} samples")

    except Exception as e:
        print(f"Error processing BioSample RunInfo: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
