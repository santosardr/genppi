#!/bin/bash
# Installation script for GenPPI Python interface

# Check if pip is available
if ! command -v pip &> /dev/null; then
    echo "Error: pip is not installed. Please install pip first."
    exit 1
fi

# Install the package in development mode
pip install -e .

echo "GenPPI Python interface installed successfully!"
echo "You can now use the 'genppi' command from the command line."