# Package requirements:
```
conda create -n bscripts python=2.7
source activate bscripts
conda install --yes numpy pandas matplotlib seaborn tqdm jupyter
conda install -c bioconda gffutils=0.8.7.1


python setup.py build
python setup.py install
```
