#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-05-10
# @DESCRIPTION: Parse JSON to get BioProject ID and GSE ID

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

import typer


def parse_json(json_file: str):
    """Parse the JSON file to get BioProject ID and GSE ID"""
    with open(json_file, "r") as f:
        data = json.load(f)

    # Extract BioProject ID and GSE ID
    document_summary = data["DocumentSummarySet"]["DocumentSummary"][0]
    # print(json.dumps(document_summary, indent=2))
    bioproject_id = document_summary.get("BioProject", "")
    gse_id = document_summary.get("Accession", "")

    return bioproject_id, gse_id


def run_command(cmd: str, verbose: bool = False):
    """Run the generated command"""
    if verbose:
        print(f"Executing command: {cmd}")

    try:
        subprocess.run(cmd, shell=True, check=True)
        print("Command executed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        sys.exit(1)


def create_output_files(
    json_file: str,
    bioproject_id: str,
    gse_id: str,
    output_dir: Optional[str] = None,
    run_cmd: bool = False,
    verbose: bool = False,
):
    """Create output files"""
    # Get directory from json_file or use provided output_dir
    dirname = output_dir if output_dir else os.path.dirname(json_file)

    # Create esearch command with full path
    output_file = os.path.join(dirname, f"{gse_id}.SraRunTable")
    cmd = f"esearch -db sra -query {bioproject_id} | efetch -format runinfo | awk 'NR==1 || /RNA-Seq/'>{output_file}"

    # Write command to shell script
    shell_script_path = os.path.join(dirname, f"00.edirect.sra.{gse_id}.sh")
    with open(shell_script_path, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write(f"# BioProject {bioproject_id}\n")
        f.write(f"# GSE {gse_id}\n\n")
        f.write(f"{cmd}\n")

    # Make shell script executable
    os.chmod(shell_script_path, 0o755)

    print(f"Created shell script: {shell_script_path}")
    print(f"When executed, the script will create: {output_file}")

    # Create biosample command
    biosample_output_file = os.path.join(
        dirname, f"{gse_id}.edirect.biosample.runinfo"
    )
    biosample_cmd = f"esearch -db sra -query {bioproject_id} | elink -target biosample | efetch -format runinfo >{biosample_output_file}"

    # Write biosample command to shell script
    biosample_shell_script_path = os.path.join(
        dirname, f"00.edirect.biosample.{gse_id}.sh"
    )
    with open(biosample_shell_script_path, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write(f"# BioProject {bioproject_id}\n")
        f.write(f"# GSE {gse_id}\n\n")
        f.write(f"{biosample_cmd}\n")

    # Make biosample shell script executable
    os.chmod(biosample_shell_script_path, 0o755)

    print(f"Created biosample shell script: {biosample_shell_script_path}")
    print(f"When executed, the script will create: {biosample_output_file}")

    # Run the commands if requested
    if run_cmd:
        print("Running commands...")
        run_command(cmd, verbose)
        run_command(biosample_cmd, verbose)


app = typer.Typer(
    help="Parse JSON to get BioProject ID and GSE ID and create SRA run scripts"
)


@app.command()
def main(
    json_file: Path = typer.Argument(
        ...,
        help="Path to the input JSON file containing BioProject and GSE information",
    ),
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Output directory for SraRunTable and shell script (defaults to same directory as JSON)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Print verbose information during processing",
    ),
    run: bool = typer.Option(
        False,
        "--run",
        "-r",
        help="Execute the generated command after creating the shell script",
    ),
):
    """
    Parse JSON to get BioProject ID and GSE ID and create SRA run scripts.
    The script generates shell files to fetch SRA run table and BioSample information.
    """
    # Check if the file exists
    if not json_file.exists():
        print(f"Error: File {json_file} does not exist")
        sys.exit(1)

    # Check if output directory exists and create if specified
    if output_dir and not output_dir.exists():
        try:
            output_dir.mkdir(parents=True)
            if verbose:
                print(f"Created output directory: {output_dir}")
        except OSError:
            print(f"Error: Could not create output directory: {output_dir}")
            sys.exit(1)

    # Parse JSON
    if verbose:
        print(f"Parsing JSON file: {json_file}")
    bioproject_id, gse_id = parse_json(str(json_file))

    if not bioproject_id or not gse_id:
        print(
            f"Error: Could not extract BioProject ID or GSE ID from {json_file}"
        )
        sys.exit(1)

    if verbose:
        print(f"Found BioProject ID: {bioproject_id}")
        print(f"Found GSE ID: {gse_id}")

    # Create output files and optionally run the command
    create_output_files(
        str(json_file),
        bioproject_id,
        gse_id,
        str(output_dir) if output_dir else None,
        run,
        verbose,
    )


if __name__ == "__main__":
    app()
