from biopandas.pdb import PandasPdb
import numpy as np
import json

PATH_XLS_INTRA_NONRED = 'XLs_intra_nonred.csv'
PATH_WATERSHED_LABELS = 'watershed_labels.json'

UNRELAXED_PDB_SUFFIX = '-unrelaxed_model_1.pdb'
RELAXED_PDB_SUFFIX = '-relaxed_model_1.pdb'
PTM_PDB_SUFFIX = '-relaxed_model_1_ptm.pdb'
RESULTS_PKL_SUFFIX = '-result_model_1_ptm.pkl'
PAE_PKL_SUFFIX = '-pae_model_1_ptm.pkl'

# Lists of hex colors for labels

colors9 = [
    '#440154',  # purple
    '#482677',  # indigo
    '#39568C',  # blue
    '#2D708E',  # blue-teal
    '#238A8D',  # teal
    '#29AF7F',  # green-teal
    '#73D055',  # green
    '#B8DE29',  # yellow-green
    '#FDE725',  # yellow
]

colors1 = [colors9[0]]
colors2 = [colors9[0], colors9[8]]
colors3 = [colors9[0], colors9[4], colors9[8]]
colors4 = [colors9[0], colors9[2], colors9[5], colors9[8]]
colors5 = [colors9[0], colors9[2], colors9[4], colors9[6], colors9[8]]
colors6 = [colors9[0], colors9[1], colors9[3], colors9[4], colors9[6], colors9[8]]
colors7 = colors9[0:4] + colors9[5:8]
colors8 = colors9[0:4] + colors9[5:9]

color_palettes = {
    1: colors1,
    2: colors2,
    3: colors3,
    4: colors4,
    5: colors5,
    6: colors6,
    7: colors7,
    8: colors8,
    9: colors9
}


pdbs = {}


def sanitize_protein(protein_id):
    """Replace / in protein ID with _ for compatibility with file system"""
    return protein_id.replace("/", "_")


def calculate_dist(pdb, xl):
    """Calculate the distance between residues in the PDB based on the crosslink"""
    # Find the C atom for each residue
    res_1 = pdb.loc[pdb['residue_number'] == xl.Res1]
    res_2 = pdb.loc[pdb['residue_number'] == xl.Res2]

    res_1_coords = res_1[['x_coord', 'y_coord', 'z_coord']].to_numpy()
    res_2_coords = res_2[['x_coord', 'y_coord', 'z_coord']].to_numpy()

    return np.linalg.norm(res_1_coords - res_2_coords)


def read_pdb_ca(pdb_filepath):
    ppdb = PandasPdb().read_pdb(pdb_filepath)  # reads pdb into pandas df
    pdbatom = ppdb.df['ATOM']  # extracts only the atom data type from pdb
    return pdbatom[(pdbatom['atom_name'] == 'CA')]  # selects the only the CA atoms


def read_protein_pdb(protein, pdb_suffix):
    if pdb_suffix not in pdbs:
        pdbs[pdb_suffix] = {}
    elif protein in pdbs[pdb_suffix]:
        return pdbs[pdb_suffix][protein]
    pdbs_with_suffix = pdbs[pdb_suffix]

    pdb_filepath = 'models/' + sanitize_protein(protein) + pdb_suffix
    pdb_ca = read_pdb_ca(pdb_filepath)
    if pdb_ca is None:
        error('PDB not found: ' + pdb_filepath)

    pdbs_with_suffix[protein] = pdb_ca
    return pdb_ca


def fill_xls_from_pdb(protein_name, threshold, pdb, xls):
    xls['af_distance'] = xls.apply(
        lambda x: x.af_distance if x.Protein1 != protein_name else calculate_dist(pdb, x), axis=1)
    return xls[(xls['af_distance'] < threshold) & (xls['Protein1'] == protein_name)]



import pandas as pd
from chimerax.core.commands import run

CROSSLINK_THRESHOLD = 30


def split_xl_agreement_violation_pairs(protein_name, xls):
    protein_xls = xls.loc[xls['Protein1'] == protein_name]
    protein_xls_agreed = protein_xls[protein_xls['af_distance'] <= CROSSLINK_THRESHOLD]
    protein_xls_violated = protein_xls[protein_xls['af_distance'] > CROSSLINK_THRESHOLD]
    xl_positions_agreed = protein_xls_agreed[['newRes1', 'newRes2']]
    xl_positions_violated = protein_xls_violated[['newRes1', 'newRes2']]
    xl_agreements = [tuple(x) for x in xl_positions_agreed.to_numpy()]
    xl_violations = [tuple(x) for x in xl_positions_violated.to_numpy()]
    return xl_agreements, xl_violations


XLs_intra_nonred = pd.read_csv(PATH_XLS_INTRA_NONRED)
XLs_intra_nonred.loc[:, 'af_distance'] = None

with open(PATH_WATERSHED_LABELS) as watershed_labels_file:
    labels = json.load(watershed_labels_file)

proteins = ['IFT172', 'IDP2\xa0', 'Q23TY1', 'Q234E6', 'I7MHP2', 'CDC27', 'SARS', 'Q24GN6', 'DLD', 'EEF2',
                'Q24FH6', 'GNPDA1/2', 'Q23A15', 'PGK1', 'SPEF2', 'EIF4A', 'Q22F29', 'Q23DV1', 'Q22T19', 'I7LWT9',
                'I7ML23', 'ILS1', 'I7LUZ1', 'DARS', 'Q24CJ0', 'CASC1', 'BBC52', 'CAPN-1', 'I7M2E8', 'I7M2A9', 'EPC1',
                'KARS', 'I7LXF9', 'I7LX35', 'GRS1', 'I7MCM4', 'EARS', 'YARS', 'I7MJ59', 'I7MKN3', 'TFA', 'I7M3K6',
                'LRS1', 'PRKAR1A', 'RANGAD1', 'BBC118', 'HSPA4', 'IGR3', 'MRNC57', 'VMA1', 'I7MHD4', 'SPAG17',
                'DPY30-1', 'MRNO36', 'RPS0', 'I7MFS4', 'Q22MP6', 'jacalin-1', 'Q22AS9', 'MPK3', 'DIC3', 'Q24C62',
                'DIC2', 'I7M6H8', 'EPA1', 'PGM1', 'UBXN2A/B', 'FTT18', 'CCDC81-1', 'PGI1', 'EMAP5/6-2', 'VMA5',
                'EMAP5/6-1', 'CCDC96', 'TTC18', 'DYH24', 'I7LZI8', 'NPEPL1/LAP3', 'I7MEJ5', 'EF1A', 'PRKG1', 'BBC73',
                'ACBD7', 'RACK1', 'I7M328', 'TRAF3IP1', 'Q23F83', 'I7M350', 'AK1-1', 'I7MDK2', 'I7M0R3', 'HSP90AA1',
                'FBPA', 'ACO1-2', 'ME1/2/3', 'I7LW80', 'TRS1']


for protein in proteins:
    # Get the XL violations
    try:
        pdb = read_protein_pdb(protein, PTM_PDB_SUFFIX)
    except FileNotFoundError:
        continue
    fill_xls_from_pdb(protein, CROSSLINK_THRESHOLD, pdb, XLs_intra_nonred)

    agreed_xl_pairs, violated_xl_pairs = split_xl_agreement_violation_pairs(protein, XLs_intra_nonred)

    # Chimera X
    run(session, 'close session')
    run(session, 'log clear')

    ptm_pdb_filepath = 'models/' + sanitize_protein(protein) + PTM_PDB_SUFFIX
    try:
        protein_model = run(session, 'open ' + ptm_pdb_filepath)
    except FileNotFoundError:
        continue

    unrelaxed_pdb_filepath = 'models/' + sanitize_protein(protein) + UNRELAXED_PDB_SUFFIX
    opened_models = run(session, 'open ' + unrelaxed_pdb_filepath)

    run(session, 'matchmaker #2@ca to #1@ca')
    run(session, 'hide #2 models')

    # Color by labels
    run(session, 'color white')
    run(session, 'transparency #2 75 All')

    protein_labels = labels[sanitize_protein(protein)]
    unique_labels, unique_indices, unique_counts = np.unique(protein_labels, return_index=True, return_counts=True)
    num_labels = np.count_nonzero(unique_labels)
    color_palette = color_palettes.get(num_labels)

    for i in range(len(unique_labels)):
        label = unique_labels[i]
        if label == 0:
            continue
        color = color_palette[label - 1]
        start = unique_indices[i]
        end = start + unique_counts[i]
        run(session, 'color #1,2:{}-{} {}'.format(start, end, color))

    # Draw violations
    for pair in violated_xl_pairs:
        if protein_labels[pair[0]] == protein_labels[pair[1]] and protein_labels[pair[0]] != 0:
            color = 'magenta'
        else:
            color = 'darkred'
        run(session, 'distance #1:{}@ca #1:{}@ca dashes 0 radius 0.5 color {}'.format(pair[0], pair[1], color))
    # run(session, 'distance style dashes 0 radius 0.2 color black')

    # Draw agreements
    for pair in agreed_xl_pairs:
        if protein_labels[pair[0]] == protein_labels[pair[1]] and protein_labels[pair[0]] != 0:
            color = 'cyan'
        else:
            color = 'darkblue'
        run(session, 'distance #1:{}@ca #1:{}@ca dashes 0 radius 0.5 color {}'.format(pair[0], pair[1], color))

    run(session, 'set bgColor white')
    run(session, 'lighting soft depthCue false')
    run(session, 'lighting model #3 directional false shadows false')
    run(session, 'view clip false')
    run(session, 'hide #3.1 models')

    # Save Results
    run(session, 'save chx/PTM/Images/{}.png supersample 3'.format(sanitize_protein(protein)))
    run(session, 'save chx/PTM/{}.cxs format session'.format(sanitize_protein(protein)))
    run(session, 'log save chx/PTM/{}_log.html'.format(sanitize_protein(protein)))
