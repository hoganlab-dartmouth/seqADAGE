# seqADAGE

> ADAGE models trained on RNAseq compendia of *Pseduomonas aeruginosa* gene expression.

## Setup

This project requires the keras library for python (python 3 as of 6/22/21) and includes a couple class definitions designed for model inspection and interpretation. 

I've made an up-to-date environment for this project (tfk.yml)

name: tfk
channels:
  - defaults
  - anaconda
dependencies:
  - jupyter
  - matplotlib
  - numpy
  - python=3.7
  - scikit-learn
  - tensorflow
  - keras
  - pandas
  - matplotlib
  - seaborn
  - numpy
  - scikit-learn
  - scipy
  - statsmodels

## How to train a seqADAGE model

See jupyter notebook for test case example of training a model using run_model.py and assessing it breifly using pathway_enrichment.R (ref adage).

## Documentation

Vignettes and examples to come.

## FAQ

Information pointing to compendium creation to come.

## Support

## License 

