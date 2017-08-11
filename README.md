[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

# Iudex 
Iudex is an analysis pipeline for the mapping of genetrap insertion events in haploid human cells using Illumina next generation sequencing. See the methods outlined in [Davis *et. al.*](http://www.cell.com/cell-reports/fulltext/S2211-1247(15)00552-5?_returnURL=http%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124715005525%3Fshowall%3Dtrue), or any of the high throughput screening experiments from the [Brummelkamp Lab](https://www.nki.nl/divisions/biochemistry/brummelkamp-t-group/publications-brummelkamp/).

Libraires sequenced using these gene trap plasmids and with linear PCR primers are designed such that the first base in the read will mark the retroviral insertion site. This is what *Iudex* will use to obtain the count data for insertion events in a given gene thus annotating gene hits.


More Information
----------------
  - [Software Requirements](https://github.com/davisem/Iudex/blob/master/docs/software_requirements.md)
  - [Installation](https://github.com/davisem/Iudex/blob/master/docs/installation.md)
  - [Usage](https://github.com/davisem/Iudex/blob/master/docs/usage.md)
  - [Performance Tuning](https://github.com/davisem/Iudex/blob/master/docs/perf.md)