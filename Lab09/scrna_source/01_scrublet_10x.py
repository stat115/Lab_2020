import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import gzip
import csv

x10_dr = sys.argv[1]
print("Parsing 10x scRNA data from: " + x10_dr)


input_dir = x10_dr + "/"

# import cells
cells=[]
f = gzip.open(input_dir + "barcodes.tsv.gz", 'rt')
for line in f.readlines():
    cells.append(line.strip())
f.close()

# Import counts and run scrublet
counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx.gz').T.tocsc()
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
rows = zip(cells, doublet_scores, predicted_doublets)
out_file_path = "scrublet_out/" + x10_dr + ".scrub.tsv"
with open(out_file_path, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["barcode", "score", "called"])

    for row in rows:
        writer.writerow(row)