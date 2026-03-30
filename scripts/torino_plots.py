import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

def plot_isoform_annotations(annotation_exons, gene, colores=None, start=None, end=None, figsize=None, lwidth=5, iso_order=None, axes=None, xlim=None):

    gene_exons = annotation_exons.loc[annotation_exons.gene_id == gene]

    if iso_order is None:
        isoforms = sorted(gene_exons.transcript_id.unique())

    else:
        isoforms = iso_order
    # else:
    #     isoforms = [gene + '.' + x for x in iso_order]

    isoform_dict = {}
    for i, iso in enumerate(isoforms):
        isoform_name = f'isoform_{str(i+1)}'
        df = gene_exons.loc[gene_exons.transcript_id == iso].copy()
        df['transcript_id'] =  f'{gene}.{isoform_name}'
        isoform_dict.update({isoform_name:{'df':df}})

    chrom = list(annotation_exons.chrom)[0]
    if start is None:
        start = str(np.min([int(list(gene_exons.start)[0]), int(list(gene_exons.start)[0])]) - 1000)
    if end is None:
        end = str(np.max([int(list(gene_exons.end)[-1]), int(list(gene_exons.end)[-1])]) + 1000)

    coords = [f'{chrom}:{start}', f'{chrom}:{end}']

    # print(isoform_dict)

    plot_gene_isoforms(isoform_dict, coords, color_list = colores, figsize=figsize, lwidth=lwidth, axes=axes, xlim=xlim)


def plot_gene_isoforms(isoforms_dict, coordinates, color_list = None, axes=None, figsize=None, lwidth=5, xlim=None):

    if xlim is None:
        xlim1 = int(coordinates[0].split(':')[1])
        xlim2 = int(coordinates[-1].split(':')[1])
        xlim = [xlim1, xlim2]
    else:
        xlim1 = xlim[0]
        xlim2 = xlim[1]
        xlim = (xlim1, xlim2)
    

    if color_list is None:
        color_list = sns.color_palette("tab10")

    K = len(isoforms_dict)

    if figsize is None:
        figsize=(20, 3)

    if axes is None:
        fig, axes = plt.subplots(K, 1, figsize=figsize)

    for i in range(K):
        isoform_df = isoforms_dict[f'isoform_{str(i+1)}']['df']
        if K == 1:
            ax = axes
        else:
            ax = axes[i]
        color = color_list[i]
        plot_isoform(isoform_df, ax, color, xlim[0], xlim[1], lwidth=lwidth)


def plot_isoform(isoform_df, ax, color, xlim1, xlim2, lwidth):
    # print(isoform_df)
    is_first = True
    for idx, row in isoform_df.iterrows():
        start = int(row.start)
        end = int(row.end)
        if is_first:
            first = end
            is_first = False
    
        ax.fill_between([start, end], [0, 0], [1, 1], color = color, zorder=2)
    
    ax.plot([first, start], [0.5, 0.5], c=color, linewidth=lwidth)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[['bottom', 'top', 'right', 'left']].set_visible(False)
    ax.set_xlim([xlim1, xlim2])

    # print(xlim1, xlim2)

gtex_colors = {
  "Adipose - Subcutaneous": {
    "tissue_abbrv": "ADPSBQ", 
    "tissue_color_hex": "FFA54F", 
    "tissue_color_rgb": "255,165,79"
  }, 
  "Adipose - Visceral (Omentum)": {
    "tissue_abbrv": "ADPVSC", 
    "tissue_color_hex": "EE9A00", 
    "tissue_color_rgb": "238,154,0"
  }, 
  "Adrenal Gland": {
    "tissue_abbrv": "ADRNLG", 
    "tissue_color_hex": "8FBC8F", 
    "tissue_color_rgb": "143,188,143"
  }, 
  "Artery - Aorta": {
    "tissue_abbrv": "ARTAORT", 
    "tissue_color_hex": "8B1C62", 
    "tissue_color_rgb": "139,28,98"
  }, 
  "Artery - Coronary": {
    "tissue_abbrv": "ARTCRN", 
    "tissue_color_hex": "EE6A50", 
    "tissue_color_rgb": "238,106,80"
  }, 
  "Artery - Femoral": {
    "tissue_abbrv": "ARTFMR", 
    "tissue_color_hex": "FF4500", 
    "tissue_color_rgb": "255,69,0"
  }, 
  "Artery - Tibial": {
    "tissue_abbrv": "ARTTBL", 
    "tissue_color_hex": "FF0000", 
    "tissue_color_rgb": "255,0,0"
  }, 
  "Bladder": {
    "tissue_abbrv": "BLDDER", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Brain - Amygdala": {
    "tissue_abbrv": "BRNAMY", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Anterior cingulate cortex (BA24)": {
    "tissue_abbrv": "BRNACC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Caudate (basal ganglia)": {
    "tissue_abbrv": "BRNCDT", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Cerebellar Hemisphere": {
    "tissue_abbrv": "BRNCHB", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Cerebellum": {
    "tissue_abbrv": "BRNCHA", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Cortex": {
    "tissue_abbrv": "BRNCTXA", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Frontal Cortex (BA9)": {
    "tissue_abbrv": "BRNCTXB", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Hippocampus": {
    "tissue_abbrv": "BRNHPP", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Hypothalamus": {
    "tissue_abbrv": "BRNHPT", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Nucleus accumbens (basal ganglia)": {
    "tissue_abbrv": "BRNNCC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Putamen (basal ganglia)": {
    "tissue_abbrv": "BRNPTM", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Spinal cord (cervical c-1)": {
    "tissue_abbrv": "BRNSPC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Substantia nigra": {
    "tissue_abbrv": "BRNSNG", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Breast - Mammary Tissue": {
    "tissue_abbrv": "BREAST", 
    "tissue_color_hex": "00CDCD", 
    "tissue_color_rgb": "0,205,205"
  }, 
  "Cells - EBV-transformed lymphocytes": {
    "tissue_abbrv": "LCL", 
    "tissue_color_hex": "EE82EE", 
    "tissue_color_rgb": "238,130,238"
  }, 
  "Cells - Transformed fibroblasts": {
    "tissue_abbrv": "FIBRBLS", 
    "tissue_color_hex": "9AC0CD", 
    "tissue_color_rgb": "154,192,205"
  }, 
  "Cervix - Ectocervix": {
    "tissue_abbrv": "CVXECT", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Cervix - Endocervix": {
    "tissue_abbrv": "CVSEND", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Colon - Sigmoid": {
    "tissue_abbrv": "CLNSGM", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Colon - Transverse": {
    "tissue_abbrv": "CLNTRN", 
    "tissue_color_hex": "EEC591", 
    "tissue_color_rgb": "238,197,145"
  }, 
  "Esophagus - Gastroesophageal Junction": {
    "tissue_abbrv": "ESPGEJ", 
    "tissue_color_hex": "8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus - Mucosa": {
    "tissue_abbrv": "ESPMCS", 
    "tissue_color_hex": "8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus - Muscularis": {
    "tissue_abbrv": "ESPMSL", 
    "tissue_color_hex": "CDAA7D", 
    "tissue_color_rgb": "205,170,125"
  }, 
  "Fallopian Tube": {
    "tissue_abbrv": "FLLPNT", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Heart - Atrial Appendage": {
    "tissue_abbrv": "HRTAA", 
    "tissue_color_hex": "B452CD", 
    "tissue_color_rgb": "180,82,205"
  }, 
  "Heart - Left Ventricle": {
    "tissue_abbrv": "HRTLV", 
    "tissue_color_hex": "7A378B", 
    "tissue_color_rgb": "122,55,139"
  }, 
  "Kidney - Cortex": {
    "tissue_abbrv": "KDNCTX", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Kidney - Medulla": {
    "tissue_abbrv": "KDNMDL", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Liver": {
    "tissue_abbrv": "LIVER", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Lung": {
    "tissue_abbrv": "LUNG", 
    "tissue_color_hex": "9ACD32", 
    "tissue_color_rgb": "154,205,50"
  }, 
  "Minor Salivary Gland": {
    "tissue_abbrv": "SLVRYG", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Muscle - Skeletal": {
    "tissue_abbrv": "MSCLSK", 
    "tissue_color_hex": "7A67EE", 
    "tissue_color_rgb": "122,103,238"
  }, 
  "Nerve - Tibial": {
    "tissue_abbrv": "NERVET", 
    "tissue_color_hex": "FFD700", 
    "tissue_color_rgb": "255,215,0"
  }, 
  "Ovary": {
    "tissue_abbrv": "OVARY", 
    "tissue_color_hex": "FFB6C1", 
    "tissue_color_rgb": "255,182,193"
  }, 
  "Pancreas": {
    "tissue_abbrv": "PNCREAS", 
    "tissue_color_hex": "CD9B1D", 
    "tissue_color_rgb": "205,155,29"
  }, 
  "Pituitary": {
    "tissue_abbrv": "PTTARY", 
    "tissue_color_hex": "B4EEB4", 
    "tissue_color_rgb": "180,238,180"
  }, 
  "Prostate": {
    "tissue_abbrv": "PRSTTE", 
    "tissue_color_hex": "D9D9D9", 
    "tissue_color_rgb": "217,217,217"
  }, 
  "Skin - Not Sun Exposed (Suprapubic)": {
    "tissue_abbrv": "SKINNS", 
    "tissue_color_hex": "3A5FCD", 
    "tissue_color_rgb": "58,95,205"
  }, 
  "Skin - Sun Exposed (Lower leg)": {
    "tissue_abbrv": "SKINS", 
    "tissue_color_hex": "1E90FF", 
    "tissue_color_rgb": "30,144,255"
  }, 
  "Small Intestine - Terminal Ileum": {
    "tissue_abbrv": "SNTTRM", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Spleen": {
    "tissue_abbrv": "SPLEEN", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Stomach": {
    "tissue_abbrv": "STMACH", 
    "tissue_color_hex": "FFD39B", 
    "tissue_color_rgb": "255,211,155"
  }, 
  "Testis": {
    "tissue_abbrv": "TESTIS", 
    "tissue_color_hex": "A6A6A6", 
    "tissue_color_rgb": "166,166,166"
  }, 
  "Thyroid": {
    "tissue_abbrv": "THYROID", 
    "tissue_color_hex": "008B45", 
    "tissue_color_rgb": "0,139,69"
  }, 
  "Uterus": {
    "tissue_abbrv": "UTERUS", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Vagina": {
    "tissue_abbrv": "VAGINA", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Whole Blood": {
    "tissue_abbrv": "WHLBLD", 
    "tissue_color_hex": "FF00FF", 
    "tissue_color_rgb": "255,0,255"
  }
}




gtex_colors_cereb = {
  "Adipose - Subcutaneous": {
    "tissue_abbrv": "ADPSBQ", 
    "tissue_color_hex": "FFA54F", 
    "tissue_color_rgb": "255,165,79"
  }, 
  "Adipose - Visceral (Omentum)": {
    "tissue_abbrv": "ADPVSC", 
    "tissue_color_hex": "EE9A00", 
    "tissue_color_rgb": "238,154,0"
  }, 
  "Adrenal Gland": {
    "tissue_abbrv": "ADRNLG", 
    "tissue_color_hex": "8FBC8F", 
    "tissue_color_rgb": "143,188,143"
  }, 
  "Artery - Aorta": {
    "tissue_abbrv": "ARTAORT", 
    "tissue_color_hex": "8B1C62", 
    "tissue_color_rgb": "139,28,98"
  }, 
  "Artery - Coronary": {
    "tissue_abbrv": "ARTCRN", 
    "tissue_color_hex": "EE6A50", 
    "tissue_color_rgb": "238,106,80"
  }, 
  "Artery - Femoral": {
    "tissue_abbrv": "ARTFMR", 
    "tissue_color_hex": "FF4500", 
    "tissue_color_rgb": "255,69,0"
  }, 
  "Artery - Tibial": {
    "tissue_abbrv": "ARTTBL", 
    "tissue_color_hex": "FF0000", 
    "tissue_color_rgb": "255,0,0"
  }, 
  "Bladder": {
    "tissue_abbrv": "BLDDER", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Brain - Amygdala": {
    "tissue_abbrv": "BRNAMY", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Anterior cingulate cortex (BA24)": {
    "tissue_abbrv": "BRNACC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Caudate (basal ganglia)": {
    "tissue_abbrv": "BRNCDT", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Cerebellar Hemisphere": {
    "tissue_abbrv": "BRNCHB", 
    "tissue_color_hex": "C4C400", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Cerebellum": {
    "tissue_abbrv": "BRNCHA", 
    "tissue_color_hex": "C4C400", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Cortex": {
    "tissue_abbrv": "BRNCTXA", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Frontal Cortex (BA9)": {
    "tissue_abbrv": "BRNCTXB", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Hippocampus": {
    "tissue_abbrv": "BRNHPP", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Hypothalamus": {
    "tissue_abbrv": "BRNHPT", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Nucleus accumbens (basal ganglia)": {
    "tissue_abbrv": "BRNNCC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Putamen (basal ganglia)": {
    "tissue_abbrv": "BRNPTM", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Spinal cord (cervical c-1)": {
    "tissue_abbrv": "BRNSPC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain - Substantia nigra": {
    "tissue_abbrv": "BRNSNG", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Breast - Mammary Tissue": {
    "tissue_abbrv": "BREAST", 
    "tissue_color_hex": "00CDCD", 
    "tissue_color_rgb": "0,205,205"
  }, 
  "Cells - EBV-transformed lymphocytes": {
    "tissue_abbrv": "LCL", 
    "tissue_color_hex": "EE82EE", 
    "tissue_color_rgb": "238,130,238"
  }, 
  "Cells - Transformed fibroblasts": {
    "tissue_abbrv": "FIBRBLS", 
    "tissue_color_hex": "9AC0CD", 
    "tissue_color_rgb": "154,192,205"
  }, 
  "Cervix - Ectocervix": {
    "tissue_abbrv": "CVXECT", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Cervix - Endocervix": {
    "tissue_abbrv": "CVSEND", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Colon - Sigmoid": {
    "tissue_abbrv": "CLNSGM", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Colon - Transverse": {
    "tissue_abbrv": "CLNTRN", 
    "tissue_color_hex": "EEC591", 
    "tissue_color_rgb": "238,197,145"
  }, 
  "Esophagus - Gastroesophageal Junction": {
    "tissue_abbrv": "ESPGEJ", 
    "tissue_color_hex": "8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus - Mucosa": {
    "tissue_abbrv": "ESPMCS", 
    "tissue_color_hex": "8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus - Muscularis": {
    "tissue_abbrv": "ESPMSL", 
    "tissue_color_hex": "CDAA7D", 
    "tissue_color_rgb": "205,170,125"
  }, 
  "Fallopian Tube": {
    "tissue_abbrv": "FLLPNT", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Heart - Atrial Appendage": {
    "tissue_abbrv": "HRTAA", 
    "tissue_color_hex": "B452CD", 
    "tissue_color_rgb": "180,82,205"
  }, 
  "Heart - Left Ventricle": {
    "tissue_abbrv": "HRTLV", 
    "tissue_color_hex": "7A378B", 
    "tissue_color_rgb": "122,55,139"
  }, 
  "Kidney - Cortex": {
    "tissue_abbrv": "KDNCTX", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Kidney - Medulla": {
    "tissue_abbrv": "KDNMDL", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Liver": {
    "tissue_abbrv": "LIVER", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Lung": {
    "tissue_abbrv": "LUNG", 
    "tissue_color_hex": "9ACD32", 
    "tissue_color_rgb": "154,205,50"
  }, 
  "Minor Salivary Gland": {
    "tissue_abbrv": "SLVRYG", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Muscle - Skeletal": {
    "tissue_abbrv": "MSCLSK", 
    "tissue_color_hex": "7A67EE", 
    "tissue_color_rgb": "122,103,238"
  }, 
  "Nerve - Tibial": {
    "tissue_abbrv": "NERVET", 
    "tissue_color_hex": "FFD700", 
    "tissue_color_rgb": "255,215,0"
  }, 
  "Ovary": {
    "tissue_abbrv": "OVARY", 
    "tissue_color_hex": "FFB6C1", 
    "tissue_color_rgb": "255,182,193"
  }, 
  "Pancreas": {
    "tissue_abbrv": "PNCREAS", 
    "tissue_color_hex": "CD9B1D", 
    "tissue_color_rgb": "205,155,29"
  }, 
  "Pituitary": {
    "tissue_abbrv": "PTTARY", 
    "tissue_color_hex": "B4EEB4", 
    "tissue_color_rgb": "180,238,180"
  }, 
  "Prostate": {
    "tissue_abbrv": "PRSTTE", 
    "tissue_color_hex": "D9D9D9", 
    "tissue_color_rgb": "217,217,217"
  }, 
  "Skin - Not Sun Exposed (Suprapubic)": {
    "tissue_abbrv": "SKINNS", 
    "tissue_color_hex": "3A5FCD", 
    "tissue_color_rgb": "58,95,205"
  }, 
  "Skin - Sun Exposed (Lower leg)": {
    "tissue_abbrv": "SKINS", 
    "tissue_color_hex": "1E90FF", 
    "tissue_color_rgb": "30,144,255"
  }, 
  "Small Intestine - Terminal Ileum": {
    "tissue_abbrv": "SNTTRM", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Spleen": {
    "tissue_abbrv": "SPLEEN", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Stomach": {
    "tissue_abbrv": "STMACH", 
    "tissue_color_hex": "FFD39B", 
    "tissue_color_rgb": "255,211,155"
  }, 
  "Testis": {
    "tissue_abbrv": "TESTIS", 
    "tissue_color_hex": "A6A6A6", 
    "tissue_color_rgb": "166,166,166"
  }, 
  "Thyroid": {
    "tissue_abbrv": "THYROID", 
    "tissue_color_hex": "008B45", 
    "tissue_color_rgb": "0,139,69"
  }, 
  "Uterus": {
    "tissue_abbrv": "UTERUS", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Vagina": {
    "tissue_abbrv": "VAGINA", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Whole Blood": {
    "tissue_abbrv": "WHLBLD", 
    "tissue_color_hex": "FF00FF", 
    "tissue_color_rgb": "255,0,255"
  }
}





gtex_colors2 = {
  "Adipose_Subcutaneous": {
    "tissue_abbrv": "ADPSBQ", 
    "tissue_color_hex": "#FFA54F", 
    "tissue_color_rgb": "255,165,79"
  }, 
  "Adipose_Visceral_Omentum": {
    "tissue_abbrv": "ADPVSC", 
    "tissue_color_hex": "#EE9A00", 
    "tissue_color_rgb": "238,154,0"
  }, 
  "Adrenal_Gland": {
    "tissue_abbrv": "ADRNLG", 
    "tissue_color_hex": "#8FBC8F", 
    "tissue_color_rgb": "143,188,143"
  }, 
  "Artery_Aorta": {
    "tissue_abbrv": "ARTAORT", 
    "tissue_color_hex": "#8B1C62", 
    "tissue_color_rgb": "139,28,98"
  }, 
  "Artery_Coronary": {
    "tissue_abbrv": "ARTCRN", 
    "tissue_color_hex": "#EE6A50", 
    "tissue_color_rgb": "238,106,80"
  }, 
  "Artery_Femoral": {
    "tissue_abbrv": "ARTFMR", 
    "tissue_color_hex": "#FF4500", 
    "tissue_color_rgb": "255,69,0"
  }, 
  "Artery_Tibial": {
    "tissue_abbrv": "ARTTBL", 
    "tissue_color_hex": "#FF0000", 
    "tissue_color_rgb": "255,0,0"
  }, 
  "Bladder": {
    "tissue_abbrv": "BLDDER", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Brain_Amygdala": {
    "tissue_abbrv": "BRNAMY", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Anterior_cingulate_cortex_BA24": {
    "tissue_abbrv": "BRNACC", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Caudate_basal_ganglia": {
    "tissue_abbrv": "BRNCDT", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Cerebellar_Hemisphere": {
    "tissue_abbrv": "BRNCHB", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Cerebellum": {
    "tissue_abbrv": "BRNCHA", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Cortex": {
    "tissue_abbrv": "BRNCTXA", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Frontal_Cortex_BA9": {
    "tissue_abbrv": "BRNCTXB", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Hippocampus": {
    "tissue_abbrv": "BRNHPP", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Hypothalamus": {
    "tissue_abbrv": "BRNHPT", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Nucleus_accumbens_basal_ganglia": {
    "tissue_abbrv": "BRNNCC", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Putamen_basal_ganglia": {
    "tissue_abbrv": "BRNPTM", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Spinal_cord_cervical_c_1": {
    "tissue_abbrv": "BRNSPC", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain_Substantia_nigra": {
    "tissue_abbrv": "BRNSNG", 
    "tissue_color_hex": "#EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Breast_Mammary_Tissue": {
    "tissue_abbrv": "BREAST", 
    "tissue_color_hex": "#00CDCD", 
    "tissue_color_rgb": "0,205,205"
  }, 
  "Cells_EBV_transformed_lymphocytes": {
    "tissue_abbrv": "LCL", 
    "tissue_color_hex": "#EE82EE", 
    "tissue_color_rgb": "238,130,238"
  }, 
  "Cells_Transformed_fibroblasts": {
    "tissue_abbrv": "FIBRBLS", 
    "tissue_color_hex": "#9AC0CD", 
    "tissue_color_rgb": "154,192,205"
  }, 
  "Cervix_Ectocervix": {
    "tissue_abbrv": "CVXECT", 
    "tissue_color_hex": "#EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Cervix_Endocervix": {
    "tissue_abbrv": "CVSEND", 
    "tissue_color_hex": "#EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Colon_Sigmoid": {
    "tissue_abbrv": "CLNSGM", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Colon_Transverse": {
    "tissue_abbrv": "CLNTRN", 
    "tissue_color_hex": "#EEC591", 
    "tissue_color_rgb": "238,197,145"
  }, 
  "Esophagus_Gastroesophageal_Junction": {
    "tissue_abbrv": "ESPGEJ", 
    "tissue_color_hex": "#8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus_Mucosa": {
    "tissue_abbrv": "ESPMCS", 
    "tissue_color_hex": "#8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus_Muscularis": {
    "tissue_abbrv": "ESPMSL", 
    "tissue_color_hex": "#CDAA7D", 
    "tissue_color_rgb": "205,170,125"
  }, 
  "Fallopian_Tube": {
    "tissue_abbrv": "FLLPNT", 
    "tissue_color_hex": "#EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Heart_Atrial_Appendage": {
    "tissue_abbrv": "HRTAA", 
    "tissue_color_hex": "#B452CD", 
    "tissue_color_rgb": "180,82,205"
  }, 
  "Heart_Left_Ventricle": {
    "tissue_abbrv": "HRTLV", 
    "tissue_color_hex": "#7A378B", 
    "tissue_color_rgb": "122,55,139"
  }, 
  "Kidney_Cortex": {
    "tissue_abbrv": "KDNCTX", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Kidney_Medulla": {
    "tissue_abbrv": "KDNMDL", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Liver": {
    "tissue_abbrv": "LIVER", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Lung": {
    "tissue_abbrv": "LUNG", 
    "tissue_color_hex": "#9ACD32", 
    "tissue_color_rgb": "154,205,50"
  }, 
  "Minor_Salivary_Gland": {
    "tissue_abbrv": "SLVRYG", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Muscle_Skeletal": {
    "tissue_abbrv": "MSCLSK", 
    "tissue_color_hex": "#7A67EE", 
    "tissue_color_rgb": "122,103,238"
  }, 
  "Nerve_Tibial": {
    "tissue_abbrv": "NERVET", 
    "tissue_color_hex": "#FFD700", 
    "tissue_color_rgb": "255,215,0"
  }, 
  "Ovary": {
    "tissue_abbrv": "OVARY", 
    "tissue_color_hex": "#FFB6C1", 
    "tissue_color_rgb": "255,182,193"
  }, 
  "Pancreas": {
    "tissue_abbrv": "PNCREAS", 
    "tissue_color_hex": "#CD9B1D", 
    "tissue_color_rgb": "205,155,29"
  }, 
  "Pituitary": {
    "tissue_abbrv": "PTTARY", 
    "tissue_color_hex": "#B4EEB4", 
    "tissue_color_rgb": "180,238,180"
  }, 
  "Prostate": {
    "tissue_abbrv": "PRSTTE", 
    "tissue_color_hex": "#D9D9D9", 
    "tissue_color_rgb": "217,217,217"
  }, 
  "Skin_Not_Sun_Exposed_Suprapubic": {
    "tissue_abbrv": "SKINNS", 
    "tissue_color_hex": "#3A5FCD", 
    "tissue_color_rgb": "58,95,205"
  }, 
  "Skin_Sun_Exposed_Lower_leg": {
    "tissue_abbrv": "SKINS", 
    "tissue_color_hex": "#1E90FF", 
    "tissue_color_rgb": "30,144,255"
  }, 
  "Small_Intestine_Terminal_Ileum": {
    "tissue_abbrv": "SNTTRM", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Spleen": {
    "tissue_abbrv": "SPLEEN", 
    "tissue_color_hex": "#CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Stomach": {
    "tissue_abbrv": "STMACH", 
    "tissue_color_hex": "#FFD39B", 
    "tissue_color_rgb": "255,211,155"
  }, 
  "Testis": {
    "tissue_abbrv": "TESTIS", 
    "tissue_color_hex": "#A6A6A6", 
    "tissue_color_rgb": "166,166,166"
  }, 
  "Thyroid": {
    "tissue_abbrv": "THYROID", 
    "tissue_color_hex": "#008B45", 
    "tissue_color_rgb": "0,139,69"
  }, 
  "Uterus": {
    "tissue_abbrv": "UTERUS", 
    "tissue_color_hex": "#EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Vagina": {
    "tissue_abbrv": "VAGINA", 
    "tissue_color_hex": "#EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Whole_Blood": {
    "tissue_abbrv": "WHLBLD", 
    "tissue_color_hex": "#FF00FF", 
    "tissue_color_rgb": "255,0,255"
  }
}