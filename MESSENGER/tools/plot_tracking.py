"""Plot the current corrected tracking results."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).with_name("plot_corrected.py")),run_name="__main__")
