# QuantumESPRESSO modules wrapper

import os
import sys
import subprocess
import argparse
import importlib.util
from pathlib import Path
import stat


def qe_run():
    parser = argparse.ArgumentParser(
        description="""
                    Run QuantumESPRESSO modules in BOFS-1 with specified MOF and config. 
                    Must be run from the BOFS-1 directory.
                    """,
        epilog="""
                Usage:
                    ./qe_run.py <module> <config_name> <mof_file>
                Examples:
                    ./qe_run pwx pwx_scf_config mofs/SIWZOO_full_n2.cif
                    ./qe_run dosx dosx_config mofs/SIWZOO_full_n2.cif
                Available modules:
                    pwx      - Plane-Wave Self-Consistent Field
                    dosx     - Density of States
                    eelsx    - Electron Energy Loss Spectra
                    hpx      - Hubbard Parameters
                    lanczosx - TDDFT-Lanczos Response
                    magnonx  - Magnon Properties
                    phx      - Phonon Properties
                    projwfcx - Projected Wavefunctions
                """,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('module', help='QE module to run')
    parser.add_argument('config', help='Config name from qe_configs')
    parser.add_argument('mof_file', help='Path to MOF CIF file')
    args = parser.parse_args()

    def load_module(module_name):
        """Load the specified QE module."""
        module_path = Path.cwd() / 'bofs1' / 'qe' / f'{module_name}.py'
        spec = importlib.util.spec_from_file_location(f'{module_name}', module_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        
        return getattr(module, module_name)

    def load_config(config_name):
        """Load configuration from qe_configs."""
        config_path = Path.cwd() / 'bofs1' / 'qe' / 'qe_configs' / f'{config_name}.py'
        spec = importlib.util.spec_from_file_location(f'{config_name}', config_path)
        config = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config)
        
        return getattr(config, config_name)

    # Use provided MOF cif file
    mof_path = str(Path(args.mof_file).absolute())
    print(f"✓ Using MOF file: {mof_path}")
    # Load the QE module
    qe_function = load_module(args.module)
    print(f"✓ Loaded module: {args.module}")
    # Load config
    config = load_config(args.config)
    print(f"✓ Loaded config: {args.config}")
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
    qe_run()
