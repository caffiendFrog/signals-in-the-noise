import os

from collections import defaultdict

import anndata
import h5py
import pandas as pd
from scipy.sparse import csc_matrix


class TenXHDF5:
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Required file not found: {filename}")

        self.validate_file_format(filename)
        self.filename = filename
        self.adata = None

    def validate_file_format(self, filename) -> bool:
        with h5py.File(filename, 'r') as f:
            all_groups = list()
            f.visit(lambda x: all_groups.append(x))
            print(all_groups)
            file_structure = defaultdict(dict)
            for g in all_groups:
                g_split = g.split('/')
                sub_grp = file_structure[g_split[0]]
                if len(g_split) == 2:
                    sub_grp[g_split[1]] = []
                if len(g_split) == 3:
                    sub_grp[g_split[1]].append(g_split[2])

            if 'matrix' not in file_structure:
                raise ValueError('Expected to see "matrix" as a group in file.')

            for expected in ['barcodes', 'data', 'features', 'indices', 'indptr', 'shape']:
                if expected not in file_structure['matrix']:
                    raise ValueError(f'Expected to see "{expected}" as a subgroup of "matrix" in file.')

            for expected in ['feature_type', 'genome', 'id', 'name']:
                if expected not in file_structure['matrix']['features']:
                    raise ValueError(f'Expected to see "{expected}" as a subgroup of "matrix-features" in file.')

            return True

    def convert_data(self):
        """
        `[()]` is h5py notation to load dataset into memory (https://docs.h5py.org/en/stable/high/dataset.html)
        :return:
        """
        with h5py.File('../assets/BRCA_GSE161529_expression.h5', 'r') as f:
            matrix = f['matrix']

            # Decode barcodes
            # -- <HDF5 dataset "barcodes": shape (332168,), type "|S200">
            # -- indicates barcodes are stored as byte strings with max length of 200 bytes
            barcodes_encoded = matrix['barcodes'][()]
            barcodes = [b.decode('utf-8') for b in barcodes_encoded]

            # Load the compressed sparse column matrix (csc_matrix) components data
            # -- <HDF5 dataset "data": shape (580377724,), type "<f4">
            # -- indicates data are 32 bit little-endian float, equivalent to np.float32
            data = matrix['data'][()]
            # -- <HDF5 dataset "indices": shape (580377724,), type "<i4">
            # -- indicates data are 32 bit little-endian int, equivalent to np.int32
            indices = matrix['indices'][()]
            # -- <HDF5 dataset "indptr": shape (332169,), type "<i4">
            # -- indicates data are 32 bit little-endian int, equivalent to np.int32
            indptr = matrix['indptr'][()]
            # -- <HDF5 dataset "shape": shape (2,), type "<i8">
            # -- indicates data are 64 bit little-endian int, equivalent to np.int64
            shape = matrix['shape'][()]

            # Load and decode features
            features = matrix['features']
            # -- <HDF5 dataset "feature_type": shape (27188,), type "|S100">
            # -- indicates data are stored as byte strings with max length of 100 byes
            feature_types = [b.decode('utf-8') for b in features['feature_type'][()]]
            # -- <HDF5 dataset "genome": shape (27188,), type "|S10">
            # -- indicates data are stored as byte strings with max length of 10 byes
            feature_genomes = [b.decode('utf-8') for b in features['genome'][()]]
            # -- <HDF5 dataset "id": shape (27188,), type "|S100">
            # -- indicates data are stored as byte strings with max length of 100 bytes
            feature_ids = [b.decode('utf-8') for b in features['id'][()]]
            # -- <HDF5 dataset "name": shape (27188,), type "|S100">
            # -- indicates data are stored as byte strings with max length of 100 bytes
            feature_names = [b.decode('utf-8') for b in features['name'][()]]

        # Reconstruct the csc_matrix
        X = csc_matrix((data, indices, indptr), shape=shape)
        # convert to csr_matrix for AnnData
        X = X.T.tocsr()

        # Create the DataFrames to build the AnnData object
        # -- https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html
        obs = pd.DataFrame(index=barcodes)
        var = pd.DataFrame({
            'gene_ids': feature_ids,
            'gene_symbols': feature_names,
            'feature_types': feature_types,
            'genome': feature_genomes,
        }, index=feature_names)

        self.adata = anndata.AnnData(X=X, obs=obs, var=var)

    def write_anndata(self, output_filename):
        if not self.adata:
            return

        self.adata.write(output_filename)
