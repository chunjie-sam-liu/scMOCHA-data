#!/usr/bin/env python
#-*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-28 11:33:36
# @DESCRIPTION:
# @VERSION: v0.0.1

import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


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
                    child_dict[child_name] = [child_dict[child_name], child_content]
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
        description='Convert XML file to JSON format.'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to input XML file',
        type=validate_file_path
    )
    parser.add_argument(
        '-o', '--output',
        help='Path to output JSON file. If not specified, will use the input filename with .json extension'
    )
    parser.add_argument(
        '-p', '--pretty',
        action='store_true',
        help='Pretty-print JSON output with indentation'
    )
    parser.add_argument(
        '--indent',
        type=int,
        default=2,
        help='Indentation level for pretty-printing (default: 2)'
    )

    return parser.parse_args()


def main():
    """
    Main function to process XML file and output JSON.
    """
    args = parse_args()

    # Set input and output files
    xml_file = args.input

    if args.output:
        json_file = args.output
    else:
        json_file = str(Path(xml_file).with_suffix('.json'))

    # Ensure output directory exists
    ensure_output_dir(json_file)

    # Parse XML to JSON
    try:
        data = parse_xml_to_json(xml_file)

        # Write to JSON file
        with open(json_file, 'w') as f:
            indent = args.indent if args.pretty else None
            json.dump(data, f, indent=indent)

        print(f"Successfully parsed {xml_file} to {json_file}")

    except Exception as e:
        print(f"Error parsing XML: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()



