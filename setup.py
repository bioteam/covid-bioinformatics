from setuptools import setup

requirements = ['biopython', 'numpy', 'PyYAML', 'tmhmm.py']

setup(
    name='covid_bio',
    use_scm_version=True,
    author="Brian Osborne",
    author_email="bosborne@alum.mit.edu",
    description="Tools for COV sequence analysis",
    packages=['covid_bio'],
    scripts=['covid_bio/download_gb_by_taxid.py', 'covid_bio/feature_to_gene_and_protein.py',
             'covid_bio/seqs_to_aligns_and_hmms.py', 'covid_bio/parse_hmmsearch_files.py',
             'covid_bio/make_json.py', 'covid_bio/annotate_to_bed.py', 'covid_bio/cov_dictionary.yaml',
             'covid_bio/cov_length_variants.yaml', 'covid_bio/vars.py'],
    python_requires='>=3.7',
    setup_requires=['setuptools_scm'],
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License"
    ]
)
