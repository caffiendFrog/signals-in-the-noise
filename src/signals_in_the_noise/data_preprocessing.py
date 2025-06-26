import sys
import os

# Add parent/src to sys.path
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), "..", "src")))

from signals_in_the_noise.utilities.doublet_scoring import GhostCells
from signals_in_the_noise.utilities.tisch2_custom_h5 import TISCH2CustomH5

if __name__ == "__main__":
    original_data = TISCH2CustomH5("../../data/archive/BRCA_GSE161529_expression.h5")
    original_data.convert_data()

    ghost_cells = GhostCells(adata=original_data.adata)
    ghost_cells.score_doublet(output_filename="../data/BRCA_GSE161529_expression_doublets.h5ad")

    print("Done!")
