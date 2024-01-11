# %%
# from pydna.primer import Primer
from pydna.dseqrecord import Dseqrecord

# U_pSNR52_Fw_1 = Primer('CGTGCGAUTCTTTGAAAAGATAATGTATGA')
# TJOS_66_P2R = Primer('ACCTGCACUTAACTAATTACATGACTCGA')
# U_pSNR52_Fw_2 = Primer('AGTGCAGGUTCTTTGAAAAGATAATGTATGA')
# TJOS_65_P1R = Primer('CACGCGAUTAACTAATTACATGACTCGA')


from teemi.lims.csv_database import get_database, get_dna_from_plate_name, get_dna_from_box_name

plasmid_plates = get_database('plasmid_plates', path= 'data/06-lims/csv_database/')

backbone = [get_dna_from_plate_name('Backbone_template - p0056_(pESC-LEU-ccdB-USER) (1).fasta', 'plasmid_plates', database_path="data/06-lims/csv_database/")]
gRNA1_template = [get_dna_from_plate_name('gRNA1_template (1).fasta', 'plasmid_plates', database_path="data/06-lims/csv_database/")]
gRNA2_template = [get_dna_from_plate_name('gRNA2_template - pESC-LEU-gRNA_CroCPR-2 (1).fasta', 'plasmid_plates', database_path="data/06-lims/csv_database/")]

vector = Dseqrecord(backbone[0], circular = True)

vector.annotations

# %%

U_pSNR52_Fw_1 = get_dna_from_box_name('U_pSNR52_Fw_1', 'primer_box', database_path="data/06-lims/csv_database/")
TJOS_66_P2R = get_dna_from_box_name('TJOS_66_P2R', 'primer_box', database_path="data/06-lims/csv_database/")
U_pSNR52_Fw_2 = get_dna_from_box_name('U_pSNR52_Fw_2', 'primer_box', database_path="data/06-lims/csv_database/")
TJOS_65_P1R = get_dna_from_box_name('TJOS_65_P1R', 'primer_box', database_path="data/06-lims/csv_database/")

# %%

# Print the primers (note how the first reverse primer starts with ACCTGCAC, which is
# the reverse complement of the first 7 bases of the second forward primer, so they can be
# ligated)
print(U_pSNR52_Fw_1.seq)
print(TJOS_66_P2R.seq)
print(U_pSNR52_Fw_2.seq)
print(TJOS_65_P1R.seq)

from pydna.amplify import pcr
gRNA1_pcr_prod = pcr(U_pSNR52_Fw_1,TJOS_66_P2R, gRNA1_template)
gRNA2_pcr_prod = pcr(U_pSNR52_Fw_2,TJOS_65_P1R, gRNA2_template)

# %%
from teemi.design.cloning import USER_enzyme
from Bio.Restriction import AsiSI
# Apply user enzyme to make 5' overhangs

gRNA1_pcr_USER = USER_enzyme(gRNA1_pcr_prod)
gRNA2_pcr_USER = USER_enzyme(gRNA2_pcr_prod)

print(gRNA1_pcr_USER.seq.__repr__())
print(gRNA2_pcr_USER.seq.__repr__())

# Open the vector
_, open_vector = vector.cut(AsiSI)

print(open_vector.seq.__repr__())

# %%
from teemi.design.cloning import nicking_enzyme

vector_nicked = Dseqrecord(nicking_enzyme(open_vector))

print(vector_nicked.seq.__repr__())

# %% circularize

plasmid = (vector_nicked + gRNA1_pcr_USER + gRNA2_pcr_USER).looped()

