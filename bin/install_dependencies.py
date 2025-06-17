#!/usr/bin/env python

import subprocess
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
requirements_in = project_root / "requirements.in"
requirements_txt = project_root / "requirements.txt"

def run(cmd, desc):
    print(f"\n--- {desc} ---")
    print(f"$ {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"{desc} failed")

def main():
    # Install pip-tools if missing
    try:
        subprocess.run(["pip-compile", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError):
        run(["pip", "install", "pip-tools"], "Installing pip-tools")

    # Compile requirements.in to requirements.txt
    run(["pip-compile", str(requirements_in)], "Compiling requirements")

    # Install compiled requirements
    run(["pip", "install", "-r", str(requirements_txt)], "Installing requirements")

    print("\nDependencies updated and installed successfully.")

if __name__ == "__main__":
    main()
