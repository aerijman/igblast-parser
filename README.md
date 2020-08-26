## Simple parser of igblast results

[Igblast](https://github.com/ncbi/igblast) allows to view the matches to the germline V, D and J genes, details at rearrangement junctions and a lot more (see [documentation](https://ncbi.github.io/igblast/)).  

The output includes a lot of information in a format that can't be directly used in statistical analyses.  

This simple script parses igblast output into a csv formatted table, easy to understand and ready to include in a statistical analysis.


### Installation

pip install igblast-parser
conda install -c bioconda igblast-parser

python setup.py install

### Usage

docs

While docs are not ready:

Command line executable in unix-like systems:
```bash
	cat <igblast.output> | igblast-parser
``` 
Pipe is not obligatory as the input could be specified with the argument `--in`   

optional argument: `--out` to specify the prefix of the output csv file   
```bash
	igblast-parser --in <igblast.output> --out <parser_output>
```

In interactive python or script
```python
	f_in = open('igblast.output','r')

	# dictionary with UMI as primary key and multiple keys/UMI
	d = igblast_parse(f)

	# pandas provides very nice tables
 	import pandas
	df = pandas.DataFrame(d).T
```


