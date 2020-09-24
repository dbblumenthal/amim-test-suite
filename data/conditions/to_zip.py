import numpy as np
import pandas as pd


def to_binary(directory):

    expr = pd.read_csv(f'{directory}/expr.csv', index_col=0).transpose(copy=True)
    expr.columns = [str(column) for column in expr.columns]
    samples = list(expr.index)
    expr.index = [sample_id for sample_id in range(len(samples))]
    expr.to_pickle(f'{directory}/expression_data.pkl', protocol=3)
    phenotypes_as_df = pd.read_csv(f'{directory}/phenotype.csv', index_col=0)
    phenotypes = np.zeros(shape=len(samples), dtype=np.uint8)
    for sample_id in range(len(samples)):
        sample = samples[sample_id]
        if phenotypes_as_df.loc[sample, 'type'] == 'control':
            phenotypes[sample_id] = 0
        else:
            phenotypes[sample_id] = 1
    np.save(f'{directory}/phenotypes.npy', phenotypes)


if __name__ == '__main__':
    directories = ['GSE3790', 'GSE30219', 'GSE75214', 'GSE112680']
    for directory in directories:
        to_binary(directory)
