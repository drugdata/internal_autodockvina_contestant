import setuptools

setuptools.setup(
    name="internal_autodockvina_contestant",
    version="0.1.0",
    url="https://github.com/drugdata/custom_celpp_contestant",

    author="Jeff Wagner",
    author_email="j5wagner@ucsd.edu",

    description="Internal autodockvina contestant for D3R CELPP competition",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages(),

    install_requires=["d3r"],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
     scripts = ['internal_autodockvina_contestant/internal_autodockvina_contestant_dock.py',
                'internal_autodockvina_contestant/internal_autodockvina_contestant_ligand_prep.py', 
                'internal_autodockvina_contestant/internal_autodockvina_contestant_protein_prep.py']
)
