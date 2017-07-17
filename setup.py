from distutils.core import setup

setup(
    name='bfx',
    version='',
    packages=['pyscripts', 'crispr', 'encode', 'rnaseq',
              'clipseq', 'general', 'pathway', 'assembly'],
    package_dir={
        'crispr':'pyscripts/crispr',
        'encode':'pyscripts/encode',
        'rnaseq':'pyscripts/rnaseq',
        'clipseq':'pyscripts/clipseq',
        'general':'pyscripts/general',
        'pathway':'pyscripts/pathway',
        'assembly':'pyscripts/assembly'
    },
    include_package_data=True,
    url='',
    license='',
    author='brianyee',
    author_email='',
    description='collection of python bioinformatics scripts',
    install_requires=[
        'numpy>=1.10',
        'pandas>=0.16',
        'matplotlib>=1.5',
        'seaborn>=0.7',
    ],
)
