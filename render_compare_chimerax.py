from os.path import exists
from chimerax.core.commands import run

PATH_XLS_INTRA_NONRED = 'XLs_intra_nonred.csv'
PATH_WATERSHED_LABELS = 'watershed_labels.json'

UNRELAXED_PDB_SUFFIX = '-unrelaxed_model_1.pdb'
RELAXED_PDB_SUFFIX = '-relaxed_model_1.pdb'
RELAXED_PTM_PDB_SUFFIX = '-relaxed_model_1_ptm.pdb'
UNRELAXED_PTM_PDB_SUFFIX = '-unrelaxed_model_1_ptm.pdb'
RESULTS_PKL_SUFFIX = '-result_model_1_ptm.pkl'
PAE_PKL_SUFFIX = '-pae_model_1_ptm.pkl'

pdbs = {}


def sanitize_protein(protein_id):
    """Replace / in protein ID with _ for compatibility with file system"""
    return protein_id.replace("/", "_")


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
    # Chimera X
    run(session, 'close session')
    run(session, 'log clear')

    relaxed_pdb_filepath = 'models/' + sanitize_protein(protein) + RELAXED_PDB_SUFFIX
    if not exists(relaxed_pdb_filepath):
        continue

    protein_model = run(session, 'open ' + relaxed_pdb_filepath)

    unrelaxed_pdb_filepath = 'models/' + sanitize_protein(protein) + UNRELAXED_PDB_SUFFIX
    opened_models = run(session, 'open ' + unrelaxed_pdb_filepath)

    run(session, 'matchmaker #2@ca to #1@ca')

    # run(session, 'color #1 #66c2a5')
    # run(session, 'color #2 #fc8d62')
    run(session, 'color #1 #80227e')
    run(session, 'color #2 #70bc75')
    run(session, 'set bgColor white')
    run(session, 'lighting soft depthCue false')
    run(session, 'graphics silhouettes true')
    run(session, 'view clip false')

    # Save Results
    run(session, 'save chx/Relaxed/Images/{}.png width 1298 height 1198 supersample 3'.format(sanitize_protein(protein)))
    run(session, 'save chx/Relaxed/{}.cxs format session'.format(sanitize_protein(protein)))
    run(session, 'log save chx/Relaxed/{}_log.html'.format(sanitize_protein(protein)))
