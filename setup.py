from setuptools import setup

requirements = ['biopython', 'numpy', 'PyYAML']

setup(
    name='covid_bioinformatics',
    use_scm_version=True,
    author="Brian Osborne",
    author_email="briano@bioteam.net",
    description="Tools for COV sequence analysis",
    packages=['covid_bioinformatics'],
    python_requires='>=3.7',
    setup_requires=['setuptools_scm'],
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3"
    ]
)
