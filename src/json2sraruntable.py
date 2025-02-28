#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-05-10
# @DESCRIPTION: Parse JSON to get BioProject ID and GSE ID

import argparse
import json
import os
import subprocess
import sys


def parse_json(json_file):
    """Parse the JSON file to get BioProject ID and GSE ID"""
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Extract BioProject ID and GSE ID
    document_summary = data['DocumentSummarySet']['DocumentSummary'][0]
    # print(json.dumps(document_summary, indent=2))
    bioproject_id = document_summary.get('BioProject', '')
    gse_id = document_summary.get('Accession', '')

    return bioproject_id, gse_id

def run_command(cmd, verbose=False):
    """Run the generated command"""
    if verbose:
        print(f"Executing command: {cmd}")

    try:
        subprocess.run(cmd, shell=True, check=True)
        print("Command executed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        sys.exit(1)

def create_output_files(json_file, bioproject_id, gse_id, output_dir=None, run_cmd=False, verbose=False):
    """Create output files"""
    # Get directory from json_file or use provided output_dir
    dirname = output_dir if output_dir else os.path.dirname(json_file)

    # Create esearch command with full path
    output_file = os.path.join(dirname, f"{gse_id}.SraRunTable")
    cmd = f"esearch -db sra -query {bioproject_id} | efetch -format runinfo>{output_file}"

    # Write command to shell script
    shell_script_path = os.path.join(dirname, f"00.edirect.sra.{gse_id}.sh")
    with open(shell_script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write(f"# BioProject {bioproject_id}\n")
        f.write(f"# GSE {gse_id}\n\n")
        f.write(f"{cmd}\n")

    # Make shell script executable
    os.chmod(shell_script_path, 0o755)

    print(f"Created shell script: {shell_script_path}")
    print(f"When executed, the script will create: {output_file}")

    # Run the command if requested
    if run_cmd:
        print("Running command...")
        run_command(cmd, verbose)

def main():
    parser = argparse.ArgumentParser(
        description='Parse JSON to get BioProject ID and GSE ID and create SRA run scripts',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('json_file',
                        help='Path to the input JSON file containing BioProject and GSE information')
    parser.add_argument('-o', '--output-dir',
                        help='Output directory for SraRunTable and shell script (defaults to same directory as JSON)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print verbose information during processing')
    parser.add_argument('-r', '--run', action='store_true',
                        help='Execute the generated command after creating the shell script')

    args = parser.parse_args()

    # Check if the file exists
    if not os.path.isfile(args.json_file):
        print(f"Error: File {args.json_file} does not exist")
        sys.exit(1)

    # Check if output directory exists and create if specified
    if args.output_dir and not os.path.exists(args.output_dir):
        try:
            os.makedirs(args.output_dir)
            if args.verbose:
                print(f"Created output directory: {args.output_dir}")
        except OSError:
            print(f"Error: Could not create output directory: {args.output_dir}")
            sys.exit(1)

    # Parse JSON
    if args.verbose:
        print(f"Parsing JSON file: {args.json_file}")
    bioproject_id, gse_id = parse_json(args.json_file)

    if not bioproject_id or not gse_id:
        print(f"Error: Could not extract BioProject ID or GSE ID from {args.json_file}")
        sys.exit(1)

    if args.verbose:
        print(f"Found BioProject ID: {bioproject_id}")
        print(f"Found GSE ID: {gse_id}")

    # Create output files and optionally run the command
    create_output_files(args.json_file, bioproject_id, gse_id, args.output_dir, args.run, args.verbose)

if __name__ == "__main__":
    main()