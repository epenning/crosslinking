from biopandas.pdb import PandasPdb
import numpy as np

PATH_XLS_INTRA_NONRED = 'XLs_intra_nonred.csv'

UNRELAXED_PDB_SUFFIX = '-unrelaxed_model_1.pdb'

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


# proteins = ['IFT172', 'IDP2\xa0', 'Q23TY1', 'Q234E6', 'I7MHP2', 'CDC27', 'SARS', 'Q24GN6', 'DLD', 'EEF2', 'Q24FH6',
#             'GNPDA1/2', 'Q23A15', 'PGK1', 'SPEF2', 'EIF4A', 'Q23DV1', 'Q22T19', 'I7LWT9', 'I7ML23', 'ILS1',
#             'I7LUZ1', 'DARS', 'Q24CJ0', 'CASC1', 'BBC52', 'CAPN-1', 'I7M2E8', 'I7M2A9', 'EPC1', 'KARS', 'I7LXF9',
#             'I7LX35', 'GRS1', 'I7MCM4', 'EARS', 'YARS', 'I7MJ59', 'I7MKN3', 'TFA', 'I7M3K6', 'LRS1', 'PRKAR1A',
#             'RANGAD1', 'BBC118', 'HSPA4', 'IGR3', 'MRNC57', 'VMA1', 'I7MHD4', 'SPAG17', 'DPY30-1', 'MRNO36', 'RPS0',
#             'I7MFS4', 'Q22MP6', 'jacalin-1', 'Q22AS9', 'MPK3', 'DIC3', 'Q24C62', 'DIC2', 'I7M6H8', 'PGM1',
#             'UBXN2A/B', 'FTT18', 'CCDC81-1', 'PGI1', 'EMAP5/6-2', 'VMA5', 'EMAP5/6-1', 'CCDC96', 'TTC18',
#             'I7LZI8', 'NPEPL1/LAP3', 'I7MEJ5', 'EF1A', 'PRKG1', 'BBC73', 'ACBD7', 'RACK1', 'I7M328', 'TRAF3IP1',
#             'Q23F83', 'I7M350', 'AK1-1', 'I7MDK2', 'I7M0R3', 'HSP90AA1', 'FBPA', 'ACO1-2', 'ME1/2/3', 'I7LW80', 'TRS1',
#             'I7LZE3', 'Q23FF1', 'I7M4K0', 'I7LTB9', 'DRH29', 'PRS1']

proteins = ['IFT172', 'IDP2\xa0', 'Q23TY1', 'Q234E6', 'I7MHP2', 'CDC27', 'SARS', 'Q24GN6', 'DLD', 'EEF2', 'Q24FH6',
            'GNPDA1/2', 'Q23A15', 'PGK1', 'SPEF2', 'EIF4A', 'Q23DV1', 'Q22T19', 'I7LWT9', 'I7ML23', 'ILS1',
            'I7LUZ1', 'DARS', 'Q24CJ0', 'CASC1', 'BBC52', 'CAPN-1', 'I7M2E8', 'I7M2A9', 'EPC1', 'KARS', 'I7LXF9',
            'I7LX35', 'GRS1', 'I7MCM4', 'EARS', 'YARS', 'I7MJ59', 'I7MKN3', 'TFA', 'I7M3K6', 'LRS1', 'PRKAR1A',
            'RANGAD1', 'BBC118', 'HSPA4', 'IGR3', 'MRNC57', 'VMA1', 'I7MHD4', 'SPAG17', 'DPY30-1', 'MRNO36', 'RPS0',
            'I7MFS4', 'Q22MP6', 'jacalin-1', 'Q22AS9', 'MPK3', 'DIC3', 'Q24C62', 'DIC2', 'I7M6H8', 'PGM1',
            'UBXN2A/B', 'FTT18', 'CCDC81-1', 'PGI1', 'EMAP5/6-2', 'VMA5', 'EMAP5/6-1', 'CCDC96', 'TTC18',
            'I7LZI8', 'NPEPL1/LAP3', 'I7MEJ5', 'EF1A', 'PRKG1', 'BBC73', 'ACBD7', 'RACK1', 'I7M328', 'TRAF3IP1',
            'Q23F83', 'I7M350', 'AK1-1', 'I7MDK2', 'I7M0R3', 'HSP90AA1', 'FBPA', 'ACO1-2', 'ME1/2/3', 'I7LW80', 'TRS1',
            'I7LZE3', 'Q23FF1', 'I7M4K0', 'I7LTB9', 'DRH29', 'PRS1']

for protein in proteins:
    # Get the XL violations
    try:
        pdb = read_protein_pdb(protein, UNRELAXED_PDB_SUFFIX)
    except FileNotFoundError:
        continue
    fill_xls_from_pdb(protein, CROSSLINK_THRESHOLD, pdb, XLs_intra_nonred)

    agreed_xl_pairs, violated_xl_pairs = split_xl_agreement_violation_pairs(protein, XLs_intra_nonred)

    # Chimera X
    run(session, 'close session')
    run(session, 'log clear')
    run(session, 'cd /Users/erin/Work/crosslinking')

    unrelaxed_pdb_filepath = 'models/' + sanitize_protein(protein) + UNRELAXED_PDB_SUFFIX
    opened_models = run(session, 'open ' + unrelaxed_pdb_filepath)

    run(session, 'color bfactor palette alphafold')

    # Draw violations
    for pair in violated_xl_pairs:
        run(session, 'distance #1:{}@ca #1:{}@ca dashes 0 radius 0.5 color darkred'.format(pair[0], pair[1]))

    # Draw agreements
    for pair in agreed_xl_pairs:
        run(session, 'distance #1:{}@ca #1:{}@ca dashes 0 radius 0.5 color darkblue'.format(pair[0], pair[1]))

    run(session, 'set bgColor white')
    run(session, 'lighting soft depthCue false')
    run(session, 'view clip false')
    run(session, 'hide #2.1 models')
    # now do the agreement distances
    run(session, 'save chx/Unrelaxed/{}.cxs format session'.format(sanitize_protein(protein)))
    run(session, 'log save chx/Unrelaxed/{}_log.html'.format(sanitize_protein(protein)))

# Caitie used darkred for violation and darkblue for satisfied
