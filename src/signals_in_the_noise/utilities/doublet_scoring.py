import os

import scrublet as scr
import scanpy as sc


class GhostCells:
    def __init__(self, *, filename=None, adata=None):
        if all([
            filename is not None,
            adata is not None,
        ]):
            raise ValueError("Please specify either a filename OR the annotated data, but not both.")
        if not any([
            filename is not None,
            adata is not None,
        ]):
            raise ValueError("Please specify either a filename OR the annotated data, but not both.")
        if filename is not None:
            if not os.path.exists(filename):
                raise FileNotFoundError(f"Required file not found: {filename}")
            self.adata = sc.read_h5ad(filename)
        if adata is not None:
            self.adata = adata

    def score_doublet(self, *, output_filename=None):
        # annotate the highly variable genes
        sc.pp.highly_variable_genes(self.adata)

        # slice highly variable genes for doublet calculation
        adata_hvg = self.adata[:, self.adata.var['highly_variable']]
        expected_doublet_rate = self.calculate_doublet_rate(adata_hvg)
        print("Running Scrublet with expected doublet rate {:.3f}...".format(expected_doublet_rate))
        scrub = scr.Scrublet(adata_hvg.X, sim_doublet_ratio=2.0, expected_doublet_rate=expected_doublet_rate, random_state=43)
        doublet_scores, _ = scrub.scrub_doublets()

        # Store the doublet score
        self.adata.obs['doublet_score'] = doublet_scores
        if output_filename is not None:
            print(f"Writing annotated data to {output_filename}...")
            self.adata.write(output_filename)
        print("Scrublet run complete.")

    def calculate_doublet_rate(self, adata):
        num_observations = adata.n_obs
        return 0.0008 * num_observations

