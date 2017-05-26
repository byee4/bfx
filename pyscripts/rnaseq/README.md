# Package requirements:
```
conda create -n bscripts python=2.7
source activate bscripts
conda install --yes numpy pandas matplotlib seaborn tqdm
```
# To plot all (to see batch effect):
```
python correlation_heatmap.py -i /projects/ps-yeolab3/bay001/maps/current_annotations/se/ -o ~/heatmap.svg -m pearson -k -s 'SE.MATS.JunctionCountOnly.txt'
```
# To plot all (normalized)
```
python correlation_heatmap.py -i /projects/ps-yeolab3/bay001/maps/current_annotations/se/ -o ~/heatmap.svg -m pearson -k -s 'SE.MATS.JunctionCountOnly.txt'
```
# To plot just significant
```
python correlation_heatmap.py -i /projects/ps-yeolab3/bay001/maps/current_normed_annotations/se/ -o ~/heatmap.svg -m pearson -k -s 'SE.MATS.JunctionCountOnly.positive.txt'

python correlation_heatmap.py -i /projects/ps-yeolab3/bay001/maps/current_normed_annotations/se/ -o ~/heatmap.svg -m pearson -k -s 'SE.MATS.JunctionCountOnly.negative.txt'
```
# To plot significant nonredundant
```
python correlation_heatmap.py -i /projects/ps-yeolab3/bay001/maps/current_normed_annotations/se/ -o ~/heatmap.svg -m pearson -k -s 'SE.MATS.JunctionCountOnly.positive.nr.txt'

python correlation_heatmap.py -i /projects/ps-yeolab3/bay001/maps/current_normed_annotations/se/ -o ~/heatmap.svg -m pearson -k -s 'SE.MATS.JunctionCountOnly.negative.nr.txt'
```