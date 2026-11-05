"""Legacy L23 CLI; new experiments use ../run_l23.py and the shared engine."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neuronsim.archive import choose_cells
from neuronsim.archive import main as archive_main


def main(argv=None):
    # Preserve the earlier four-cell/profile defaults for existing commands.
    return archive_main(argv, legacy=True)


if __name__ == "__main__":
    main()
