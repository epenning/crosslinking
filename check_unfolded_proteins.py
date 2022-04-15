from os.path import exists


UNRELAXED_PDB_SUFFIX = '-unrelaxed_model_1.pdb'
RELAXED_PDB_SUFFIX = '-relaxed_model_1.pdb'
PTM_PDB_SUFFIX = '-relaxed_model_1_ptm.pdb'
RESULTS_PKL_SUFFIX = '-result_model_1_ptm.pkl'
PAE_PKL_SUFFIX = '-pae_model_1_ptm.pkl'


def sanitize_protein(protein_id):
    """Replace / in protein ID with _ for compatibility with file system"""
    return protein_id.replace("/", "_")


all_proteins = [('PRKG1', 38),
                ('EEF2', 35),
                ('MRNC57', 34),
                ('EMAP5/6-1', 30),
                ('I7MJ59', 30),
                ('FTT18', 25),
                ('BBC118', 25),
                ('EARS', 24),
                ('EMAP5/6-2', 22),
                ('Q23A15', 20),
                ('KARS', 18),
                ('PRKAR1A', 18),
                ('HSPA4', 17),
                ('CAPN-1', 17),
                ('I7MKN3', 17),
                ('SARS', 17),
                ('I7M6H8', 16),
                ('SPAG17', 15),
                ('EPA1', 15),
                ('I7M3K6', 15),
                ('GNPDA1/2', 15),
                ('UBXN2A/B', 14),
                ('IDP2\xa0', 14),
                ('PGM1', 14),
                ('BBC73', 14),
                ('ACO1-2', 13),
                ('jacalin-1', 13),
                ('I7LWT9', 13),
                ('CASC1', 13),
                ('I7M328', 13),
                ('AK1-1', 13),
                ('Q234E6', 13),
                ('HYDIN', 13),
                ('Q23F83', 12),
                ('Q23DV1', 12),
                ('TRAF3IP1', 12),
                ('TTC18', 12),
                ('IGR3', 12),
                ('I7LXF9', 11),
                ('Q22AS9', 11),
                ('I7M2E8', 11),
                ('Q24C62', 11),
                ('Q22T19', 11),
                ('Q23TY1', 11),
                ('ME1/2/3', 11),
                ('FBPA', 11),
                ('I7LXF9', 11),
                ('Q22T19', 11),
                ('Q24C62', 11),
                ('CCDC96', 11),
                ('FBPA', 11),
                ('ME1/2/3', 11),
                ('Q22AS9', 11),
                ('HSP90AA1', 11),
                ('I7MHP2', 11),
                ('Q23TY1', 11),
                ('I7M350', 11),
                ('I7M2E8', 11),
                ('IFT172', 11),
                ('DYH24', 10),
                ('PGK1', 10),
                ('I7MDK2', 10),
                ('LRS1', 10),
                ('PGI1', 10),
                ('TRS1', 10),
                ('EIF4A', 10),
                ('I7M2A9', 9),
                ('VMA5', 9),
                ('Q22MP6', 9),
                ('Q24FH6', 9),
                ('I7MCM4', 9),
                ('Q24CJ0', 9),
                ('SPEF2', 9),
                ('BBC52', 9),
                ('RPS0', 9),
                ('ILS1', 9),
                ('DPY30-1', 9),
                ('I7LUZ1', 9),
                ('MRNO36', 9),
                ('DIC3', 9),
                ('Q24GN6', 9),
                ('I7MFS4', 9),
                ('RACK1', 9),
                ('EF1A', 9),
                ('ACBD7', 9),
                ('RANGAD1', 9),
                ('NPEPL1/LAP3', 8),
                ('Q22F29', 8),
                ('EPC1', 8),
                ('DARS', 8),
                ('TFA', 8),
                ('I7M0R3', 8),
                ('DLD', 8),
                ('DIC2', 8),
                ('GRS1', 8),
                ('CDC27', 8),
                ('I7LW80', 8),
                ('VMA1', 8),
                ('I7LX35', 8),
                ('YARS', 8),
                ('I7ML23', 8),
                ('MPK3', 8),
                ('I7MEJ5', 7),
                ('I7MHD4', 7),
                ('CCDC81-1', 7),
                ('I7LZI8', 7)]

unfolded_proteins = []

for protein_XLs in all_proteins:
    protein = protein_XLs[0]
    pdb_filepath = 'models/' + sanitize_protein(protein) + '-unrelaxed_model_1.pdb'
    if not exists(pdb_filepath):
        print(protein)
        unfolded_proteins.append(protein)

print(unfolded_proteins)

#%%
