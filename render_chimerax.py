from biopandas.pdb import PandasPdb
from os.path import exists
import numpy as np

PATH_XLS_INTRA_NONRED = 'XLs_intra_nonred.csv'

UNRELAXED_PDB_SUFFIX = '-unrelaxed_model_1.pdb'
RELAXED_PDB_SUFFIX = '-relaxed_model_1.pdb'
PTM_PDB_SUFFIX = '-relaxed_model_1_ptm.pdb'
RESULTS_PKL_SUFFIX = '-result_model_1_ptm.pkl'
PAE_PKL_SUFFIX = '-pae_model_1_ptm.pkl'

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


def get_xl_violation_pairs(protein_name, xls):
    protein_xls = xls.loc[xls['Protein1'] == protein_name]
    protein_xls_violated = protein_xls[protein_xls['af_distance'] > CROSSLINK_THRESHOLD]
    xl_positions_violated = protein_xls_violated[['newRes1', 'newRes2']]
    return [tuple(x) for x in xl_positions_violated.to_numpy()]

XLs_intra_nonred = pd.read_csv(PATH_XLS_INTRA_NONRED)
XLs_intra_nonred.loc[:, 'af_distance'] = None

proteins = [
    'PRKG1',
    # 'EEF2',
    'MRNC57',
    'I7MJ59',
    'FTT18',
    'BBC118',
    'Q23A15',
    'KARS',
    'I7M6H8',
    'TRAF3IP1',
    'EIF4A',
    'SPEF2',
    'Q22MP6',
    'RPS0'
]


for protein in proteins:
    # Get the XL violations
    pdb = read_protein_pdb(protein, PTM_PDB_SUFFIX)
    fill_xls_from_pdb(protein, CROSSLINK_THRESHOLD, pdb, XLs_intra_nonred)

    violated_xl_pairs = get_xl_violation_pairs(protein, XLs_intra_nonred)

    # Chimera X
    run(session, 'close session')
    run(session, 'log clear')
    run(session, 'cd /Users/erin/Work/crosslinking')

    ptm_pdb_filepath = 'models/' + sanitize_protein(protein) + PTM_PDB_SUFFIX
    protein_model = run(session, 'open ' + ptm_pdb_filepath)

    unrelaxed_pdb_filepath = 'models/' + sanitize_protein(protein) + UNRELAXED_PDB_SUFFIX
    opened_models = run(session, 'open ' + unrelaxed_pdb_filepath)
    run(session, 'align #1@ca toAtoms #2@ca')
    run(session, 'hide #2 models')

    run(session, 'color bfactor palette alphafold')

    # Draw violations
    for pair in violated_xl_pairs:
        run(session, 'distance #1:{}@ca #1:{}@ca'.format(pair[0], pair[1]))

    run(session, 'set bgColor white')
    run(session, 'lighting flat')
    run(session, 'view clip false')
    run(session, 'hide #3.1 models')
    run(session, 'distance style dashes 0 radius 0.2 color black')
    run(session, 'save chx/{}.cxs format session'.format(sanitize_protein(protein)))
    run(session, 'log save chx/{}_log.html'.format(sanitize_protein(protein)))



