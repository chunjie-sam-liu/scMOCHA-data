#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-28 11:33:36
# @DESCRIPTION:
# @VERSION: v0.0.1

import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional

import typer


def element_to_dict(element):
    """
    Convert an XML Element to a Python dictionary.
    Handles nested elements, text content, and attributes.
    """
    result = {}

    # Add attributes if any
    if element.attrib:
        result.update(element.attrib)

    # Process child elements
    children = list(element)
    if children:
        child_dict = {}
        for child in children:
            child_name = child.tag
            child_content = element_to_dict(child)

            # Handle array-like elements (same tag name appears multiple times)
            if child_name in child_dict:
                if type(child_dict[child_name]) is list:
                    child_dict[child_name].append(child_content)
                else:
                    child_dict[child_name] = [
                        child_dict[child_name],
                        child_content,
                    ]
            else:
                child_dict[child_name] = child_content

        # Update result with processed children
        result.update(child_dict)

    # Handle text content if it exists and is not just whitespace
    if element.text is not None and element.text.strip():
        if children:
            result["text"] = element.text.strip()
        else:
            # If no children, use the text as the value
            return element.text.strip()

    return result


def parse_xml_to_json(xml_file):
    """
    Parse an XML file to JSON format.
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Start conversion from the root element
    result = {root.tag: element_to_dict(root)}

    return result


def validate_file_path(file_path: str) -> str:
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


app = typer.Typer(help="Convert XML file to JSON format.")


@app.command()
def convert(
    input_file: str = typer.Option(
        ..., "--input", "-i", help="Path to input XML file"
    ),
    output_file: Optional[str] = typer.Option(
        None,
        "--output",
        "-o",
        help="Path to output JSON file. If not specified, will use the input filename with .json extension",
    ),
    pretty: bool = typer.Option(
        False,
        "--pretty",
        "-p",
        help="Pretty-print JSON output with indentation",
    ),
    indent: int = typer.Option(
        2, "--indent", help="Indentation level for pretty-printing (default: 2)"
    ),
):
    """
    Convert XML file to JSON format.
    """
    # Validate input file
    try:
        validate_file_path(input_file)
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)

    # Set input and output files
    xml_file = input_file

    if output_file:
        json_file = output_file
    else:
        json_file = str(Path(xml_file).with_suffix(".json"))

    # Ensure output directory exists
    ensure_output_dir(json_file)

    # Parse XML to JSON
    try:
        data = parse_xml_to_json(xml_file)

        # Write to JSON file
        with open(json_file, "w") as f:
            indent_val = indent if pretty else None
            json.dump(data, f, indent=indent_val)

        typer.echo(f"Successfully parsed {xml_file} to {json_file}")

    except Exception as e:
        typer.echo(f"Error parsing XML: {e}", err=True)
        raise typer.Exit(code=1)


def main():
    """Entry point for the application"""
    app()


if __name__ == "__main__":
    main()
