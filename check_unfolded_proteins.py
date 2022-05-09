from os.path import exists

UNRELAXED_PDB_SUFFIX = '-unrelaxed_model_1.pdb'
RELAXED_PDB_SUFFIX = '-relaxed_model_1.pdb'
PTM_PDB_SUFFIX = '-relaxed_model_1_ptm.pdb'
RESULTS_PKL_SUFFIX = '-result_model_1_ptm.pkl'
PAE_PKL_SUFFIX = '-pae_model_1_ptm.pkl'


def sanitize_protein(protein_id):
    """Replace / in protein ID with _ for compatibility with file system"""
    return protein_id.replace("/", "_").replace("\xa0", "")


all_proteins = ['IFT172', 'IDP2\xa0', 'Q23TY1', 'Q234E6', 'I7MHP2', 'CDC27', 'SARS', 'Q24GN6', 'DLD', 'EEF2', 'Q24FH6',
                'GNPDA1/2', 'Q23A15', 'PGK1', 'SPEF2', 'EIF4A', 'Q22F29', 'Q23DV1', 'Q22T19', 'I7LWT9', 'I7ML23',
                'ILS1', 'I7LUZ1', 'DARS', 'Q24CJ0', 'CASC1', 'BBC52', 'CAPN-1', 'I7M2E8', 'I7M2A9', 'EPC1', 'KARS',
                'I7LXF9', 'I7LX35', 'GRS1', 'I7MCM4', 'EARS', 'YARS', 'I7MJ59', 'I7MKN3', 'TFA', 'I7M3K6', 'LRS1',
                'PRKAR1A', 'RANGAD1', 'BBC118', 'HSPA4', 'IGR3', 'MRNC57', 'VMA1', 'I7MHD4', 'SPAG17', 'DPY30-1',
                'MRNO36', 'RPS0', 'I7MFS4', 'Q22MP6', 'jacalin-1', 'Q22AS9', 'MPK3', 'DIC3', 'Q24C62', 'DIC2', 'I7M6H8',
                'EPA1', 'PGM1', 'UBXN2A/B', 'FTT18', 'CCDC81-1', 'PGI1', 'EMAP5/6-2', 'VMA5', 'EMAP5/6-1', 'CCDC96',
                'TTC18', 'DYH24', 'I7LZI8', 'NPEPL1/LAP3', 'I7MEJ5', 'EF1A', 'PRKG1', 'BBC73', 'ACBD7', 'RACK1',
                'I7M328', 'TRAF3IP1', 'Q23F83', 'I7M350', 'AK1-1', 'I7MDK2', 'I7M0R3', 'HSP90AA1', 'FBPA', 'ACO1-2',
                'ME1/2/3', 'I7LW80', 'TRS1', 'I7LZE3', 'Q23FF1', 'I7M4K0', 'I7LTB9', 'DRH29', 'PRS1']

print('All proteins:', len(all_proteins))
print('Unfolded proteins:')

unfolded_proteins = []
for protein in all_proteins:
    pdb_filepath = 'models/' + sanitize_protein(protein) + '-unrelaxed_model_1.pdb'
    if not exists(pdb_filepath):
        print(protein)
        unfolded_proteins.append(protein)

print(unfolded_proteins)

# %%
