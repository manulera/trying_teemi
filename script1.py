# %% Load sequences

from teemi.design.fetch_sequences import read_genbank_files
path = 'data/04-genetic_parts/G8H_CYP_CPR_PARTS/'
# Load 4 plasmids with promoters
pCPR_sites = read_genbank_files(path+'CPR_promoters_-_2022-06-06.gb')
# Load 10 main sequences
CPR_sites = read_genbank_files(path+'CPR_tCYC1_-_2022-06-06.gb')

# %% Convert them to Dseqrecord objects, I guess they should be circular
from pydna.dseqrecord import Dseqrecord
pCPR_sites = [Dseqrecord(seq, circular=True) for seq in pCPR_sites]
CPR_sites = [Dseqrecord(seq, circular=True) for seq in CPR_sites]

# We arrange them in a list of the things we want to combine:
list_of_seqs = [pCPR_sites, CPR_sites]

# %% Assembly
from teemi.design.combinatorial_design import DesignAssembly
# This will be used as a Kozak sequence (Ribosome recognition site)
kozak = [Dseqrecord('TCGGTC')]

CPR_combinatorial_library = DesignAssembly(list_of_seqs, list_of_pads = kozak , positions_of_pads =[1], overlap=35, target_tm = 55 )
CPR_combinatorial_library.primer_list_to_dataframe()

CPR_combinatorial_library.primer_list_to_dataframe().to_csv('output.tsv', sep='\t', index=False)


# %% Dummy example
from pydna.amplicon import Amplicon
from pydna.assembly import Assembly

part1 = [Dseqrecord('ATGGAGACAGCAACTTCTTCCCCTTTGCCCATTAAATCGAGGAGAAACAGCGAAAATTCTGGGTCTACTA', id='part1', name='part1')]
part2 = [Dseqrecord('CAGTTATACCGCATATGAACCCTTCTTTAGCAACACCGTTGACTGTGTCGACCATGGTAAATCAATCAAA', id='part2', name='part2')]
part3 = [Dseqrecord('TTCTAAAGAGTTTATGAAGTTGACCCCAGTTCGTATTAGAGATTTTGGTTCTCCTTTGAAAAACGTGTCC', id='part3', name='part3')]


pads = [Dseqrecord('GGGG')]

list_of_seqs = [part1, part2]

CPR_combinatorial_library = DesignAssembly(list_of_seqs, list_of_pads = pads, positions_of_pads =[1], overlap=35, target_tm = 55 )

unique_amplicons : list[Amplicon] = CPR_combinatorial_library.unique_amplicons

CPR_combinatorial_library.list_of_assemblies

assemble = Assembly(CPR_combinatorial_library.list_of_assemblies[0], limit = 30)
product = assemble.assemble_linear()[0]
print(product.seq)


# %% Longer assembly

pads = [Dseqrecord('GGGG'), Dseqrecord('CCCCC')]

list_of_seqs = [part1, part2, part3]

CPR_combinatorial_library = DesignAssembly(list_of_seqs, list_of_pads = pads, positions_of_pads =[1, 2], overlap=35, target_tm = 55 )

assemble = Assembly(CPR_combinatorial_library.list_of_assemblies[0], limit = 30)
product = assemble.assemble_linear()[0]
print(product.seq)
