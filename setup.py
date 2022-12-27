from setuptools import setup

requirements = ['biopython', 'numpy', 'PyYAML', 'pyTMHMM']

setup(
    name='covidbio',
    use_scm_version=True,
    author="Brian Osborne",
    author_email="bosborne@alum.mit.edu",
    description="Tools for COV sequence analysis",
    packages=['covidbio'],
    scripts=['covidbio/download_gb_by_taxid.py', 'covidbio/feature_to_gene_and_protein.py',
             'covidbio/seqs_to_aligns_and_hmms.py', 'covidbio/parse_hmmsearch_files.py',
             'covidbio/fast_uniprot_parser.py', 'covidbio/annotate_to_bed.py',
             'covidbio/utilities.py', 'covidbio/mash_utils.py', 'covidbio/swiss-to-fasta.py'],
    python_requires='>=3.7',
    setup_requires=['setuptools_scm'],
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License"
    ],
    include_package_data=True
)
