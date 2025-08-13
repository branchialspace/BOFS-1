# Run QuantumESPRESSO modules

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
                    Must be run from the same directory where bofs1_env.sh was executed.
                    """,
        epilog="""
                Usage:
                    ./qe_run.py <module> <config_name> <mof_file>
                Examples:
                    ./qe_run pwx pwx_scf mofs/SIWZOO_full_n2.cif
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
    parser.add_argument('config', help='Config name from configs.py')
    parser.add_argument('mof_file', help='Path to MOF CIF file')
    args = parser.parse_args()

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

    def load_config(config_name):
        """Load configuration from configs.py."""
        config_path = Path.cwd() / 'BOFS-1' / 'bofs1' / 'qe' / 'configs.py'
        # Load configs module
        spec = importlib.util.spec_from_file_location("configs", config_path)
        configs = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(configs)
        
        return getattr(configs, config_name)

    # Use MOF file as provided
    mof_path = str(Path(args.mof_file).absolute())
    print(f"✓ Using MOF file: {mof_path}")
    # Load the QE module
    qe_function = load_module(args.module)
    print(f"✓ Loaded module: qe_{args.module}")
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
