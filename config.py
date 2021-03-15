###
# configuration that should be handled manually and will differ between systems
###

from os.path import dirname, expanduser, splitext
from os import getenv
from subprocess import check_output

#BASE = HOME + '/zfp36l2_analysis'
BASE = dirname(__file__)
ENV = "zfp36l2"
N_THREADS = 4

###
# software
###
CONDA = check_output(['which', 'conda']).strip()
BIN = dirname(CONDA)
KALLISTO = BIN + '/kallisto'
PANDOC = BIN + '/pandoc'
SAMTOOLS = BIN + '/samtools'
OPENSSL = BIN + '/openssl'
UPDATED_PATH = 'PATH=' + BIN + ':$PATH'

R_PATH = BIN + '/Rscript'

###
# annotations
###
ENS_V = '93'
GENOME_V = 'Rnor_6.0'
SPECIES = 'rat'
SPECIES_SL = 'rattus_norvegicus'
SPECIES_LL = 'Rattus_norvegicus'

TRANSCRIPTOME_NAME = SPECIES_LL + '.' + GENOME_V
TRANSCRIPTOME_FA = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.all.fa'

# NOTE: bowtie creates files with this prefix
KALLISTO_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.kidx'

GENCODE_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode'
ENS_URL = 'ftp://ftp.ensembl.org/pub/release'

###
# functions
###
def source_r(base, fname):
    return UPDATED_PATH + ' OMP_NUM_THREADS=1 ' + R_PATH + ' --vanilla --default-packages=methods,stats,utils,graphics,grDevices -e \'setwd("{0}")\' -e \'source("{1}")\''.format(base, fname)

def source_rmd(base, file_name, output_name = None):
    if output_name is None:
        output_name = splitext(file_name)[0]
        output_name += '.html'
    return UPDATED_PATH + ' OMP_NUM_THREADS=1 ' + R_PATH + ' --vanilla --default-packages=methods,stats,utils -e \'.libPaths("~/R_library")\' -e \'setwd("{0}")\' -e \'rmarkdown::render("{1}", output_file = "{2}")\''.format(base, file_name, output_name)
#    return UPDATED_PATH + ' ' + R_ARG + ' OMP_NUM_THREADS=1 ' + R_PATH + ' --vanilla --default-packages=methods,stats,utils -e \'.libPaths("~/R_library")\' -e \'setwd("{0}")\' -e \'rmarkdown::render("{1}", output_file = "{2}")\''.format(base, file_name, output_name)

def get_sample_ids(fname):
    ret = []
    with open(fname, 'r') as fhandle:
        for line in fhandle:
            ret.append(line.strip("\n"))
    return ret
