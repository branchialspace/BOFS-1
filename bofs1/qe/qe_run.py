#!/usr/bin/env python3

"""
Run QuantumESPRESSO modules in BOFS-1.
Activates environment, runs specified QE module with MOF and config.
Must be run from the same directory where bofs1_env.sh was executed.
Usage
    ./qe_run.py <module> <mof_file> [--config <config_name>]
"""

import os
import sys
import subprocess
import argparse
import importlib.util
from pathlib import Path
import stat


def ensure_executable():
    """Make this script executable if it isn't already."""
    script_path = Path(__file__).resolve()
    current_perms = script_path.stat().st_mode
    # Check if already executable by owner
    if not (current_perms & stat.S_IXUSR):
        # Add execute permission for owner
        new_perms = current_perms | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        script_path.chmod(new_perms)
        print(f"✓ Made {script_path.name} executable")

def activate_environment():
    """Activate the bofs1_env virtual environment if not already active."""
    # Check if we're already in the correct environment
    venv_path = os.environ.get('VIRTUAL_ENV')
    if venv_path and 'bofs1_env' in venv_path:
        print(f"✓ Already in virtual environment: {venv_path}")
        return
    # Virtual environment should be in current directory
    venv_path = Path.cwd() / 'bofs1_env'
    # Re-execute this script with the activated environment
    python_exec = venv_path / 'bin' / 'python'
    # Build the command to re-run this script with activated environment
    cmd = [str(python_exec), __file__] + sys.argv[1:]
    # Set up environment variables
    new_env = os.environ.copy()
    new_env['VIRTUAL_ENV'] = str(venv_path)
    new_env['PATH'] = f"{venv_path / 'bin'}{os.pathsep}{new_env.get('PATH', '')}"
    # Remove PYTHONHOME if set (can interfere with venv)
    new_env.pop('PYTHONHOME', None)
    print(f"✓ Activating environment: {venv_path}")
    # Re-execute with the new environment
    result = subprocess.run(cmd, env=new_env)
    sys.exit(result.returncode)

def load_config(config_name=None):
    """Load configuration from configs.py."""
    config_path = Path.cwd() / 'BOFS-1' / 'bofs1' / 'qe' / 'configs.py'
    # Load configs module
    spec = importlib.util.spec_from_file_location("configs", config_path)
    configs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(configs)
    # If no config specified, return None
    if config_name is None:
        return None
    # Get the specified config
    config = getattr(configs, config_name)
    print(f"✓ Loaded config: {config_name}")
    return config

def load_module(module_name):
    """Dynamically load the specified QE module."""
    modules = {
        'pwx': 'BOFS-1/bofs1/qe/qe_pwx.py',
        'dosx': 'BOFS-1/bofs1/qe/qe_dosx.py',
        'eelsx': 'BOFS-1/bofs1/qe/qe_eelsx.py',
        'hpx': 'BOFS-1/bofs1/qe/qe_hpx.py',
        'lanczosx': 'BOFS-1/bofs1/qe/qe_lanczosx.py',
        'magnonx': 'BOFS-1/bofs1/qe/qe_magnonx.py',
        'phx': 'BOFS-1/bofs1/qe/qe_phx.py',
        'projwfcx': 'BOFS-1/bofs1/qe/qe_projwfcx.py'}
    module_path = Path.cwd() / modules[module_name]
    # Load the module
    spec = importlib.util.spec_from_file_location(f"qe_{module_name}", module_path)
    module = importlib.util.module_from_spec(spec)
    # Add the module directory to sys.path so imports work
    sys.path.insert(0, str(module_path.parent))
    spec.loader.exec_module(module)
    # Get the main function
    func_name = f"qe_{module_name}"
    return getattr(module, func_name)

def main():
    # ensure the script is executable
    ensure_executable()
    # ensure we're in the right environment
    activate_environment()
    parser = argparse.ArgumentParser(
        description='Run QuantumESPRESSO modules with specified MOF and config',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ./qe_run.py pwx mofs/SIWZOO_full_n2.cif --config pwx_scf
  ./qe_run.py dosx mofs/SIWZOO_full_n2.cif --config dosx_config

Available modules:
  pwx      - Plane-Wave Self-Consistent Field
  dosx     - Density of States
  eelsx    - Electron Energy Loss Spectra
  hpx      - Hubbard Parameters
  lanczosx - TDDFT-Lanczos Response
  magnonx  - Magnon Properties
  phx      - Phonon Properties
  projwfcx - Projected Wavefunctions
        """)
    parser.add_argument('module', help='QE module to run')
    parser.add_argument('mof_file', help='Path to MOF CIF file')
    parser.add_argument('--config', '-c', help='Config name from configs.py')
    args = parser.parse_args()
    # Use MOF file as provided
    mof_path = str(Path(args.mof_file).absolute())
    print(f"✓ Using MOF file: {mof_path}")
    # Load the QE module
    print(f"✓ Loading module: qe_{args.module}")
    qe_function = load_module(args.module)
    # Load config
    config = load_config(args.config)
    # Run the QE function
    print(f"\n{'='*60}")
    print(f"Running {args.module} calculation...")
    print(f"{'='*60}\n")
    try:
        qe_function(mof_path, config)
        print(f"\n{'='*60}")
        print(f"✓ {args.module} calculation completed")
        print(f"{'='*60}")
    except Exception as e:
        print(f"\n{'='*60}")
        print(f"✗ Error running {args.module}: {str(e)}")
        print(f"{'='*60}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
