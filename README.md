# Institute-for-System-Biology--Cancer-Genome-Cloud
## EDA of Sample Attributes
### Data Description

# Study of Sample Atrributes Used for Genomic Databases Development

The Genotype-Tissue Expression (GTEx) project aims to provide to the scientific community a resource with which to study human gene expression and regulation and its relationship to genetic variation. This project will collect and analyze multiple human tissues from donors who are also densely genotyped, to assess genetic variation within their genomes.
Correlations between genotype and tissue-specific gene expression levels will help identify regions of the genome that influence whether and how much a gene is expressed.

### importing libraries

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
%matplotlib inline
```

### importing csv file by Pandas

```python
df = pd.read_csv('3.csv')
```

### original shape of the dataframe
```python
df.shape
(14900, 64)
```
### list of original dataframe features
```python
df.columns
Index(['SUBJID', 'SAMPID', 'SMATSSCR', 'SMCENTER', 'SMPTHNTS', 'SMRIN', 'SMTS',
       'SMTSD', 'SMUBRID', 'SMTSISCH', 'SMTSPAX', 'SMNABTCH', 'SMNABTCHT',
       'SMNABTCHD', 'SMGEBTCH', 'SMGEBTCHD', 'SMGEBTCHT', 'SMAFRZE', 'SMGTC',
       'SME2MPRT', 'SMCHMPRS', 'SMNTRART', 'SMNUMGPS', 'SMMAPRT', 'SMEXNCRT',
       'SM550NRM', 'SMGNSDTC', 'SMUNMPRT', 'SM350NRM', 'SMRDLGTH', 'SMMNCPB',
       'SME1MMRT', 'SMSFLGTH', 'SMESTLBS', 'SMMPPD', 'SMNTERRT', 'SMRRNANM',
       'SMRDTTL', 'SMVQCFL', 'SMMNCV', 'SMTRSCPT', 'SMMPPDPR', 'SMCGLGTH',
       'SMGAPPCT', 'SMUNPDRD', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF', 'SMMPPDUN',
       'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 'SMMFLGTH', 'SME1ANTI',
       'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 'SMRRNART', 'SME1MPRT',
       'SMNUM5CD', 'SMDPMPRT', 'SME2PCTS'],
      dtype='object')
```
### finding unique values of all features and getting a sense of how many of them are there!
```
feature name, # 0f unique values
SUBJID 752
SAMPID 14900
SMATSSCR 5
SMCENTER 9
SMPTHNTS 7969
SMRIN 56
SMTS 31
SMTSD 51
SMUBRID 50
SMTSISCH 1654
SMTSPAX 1183
SMNABTCH 1418
SMNABTCHT 13
SMNABTCHD 712
SMGEBTCH 270
SMGEBTCHD 201
SMGEBTCHT 11
SMAFRZE 5
SMGTC 380
SME2MPRT 11850
SMCHMPRS 23
SMNTRART 11984
SMNUMGPS 1
SMMAPRT 11752
SMEXNCRT 12074
SM550NRM 1
SMGNSDTC 5691
SMUNMPRT 2
SM350NRM 1
SMRDLGTH 3
SMMNCPB 1
SME1MMRT 12097
SMSFLGTH 333
SMESTLBS 2
SMMPPD 12100
SMNTERRT 12097
SMRRNANM 12045
SMRDTTL 12100
SMVQCFL 12088
SMMNCV 1
SMTRSCPT 5747
SMMPPDPR 12102
SMCGLGTH 1
SMGAPPCT 1
SMUNPDRD 2
SMNTRNRT 12097
SMMPUNRT 11752
SMEXPEFF 12091
SMMPPDUN 12100
SME2MMRT 12097
SME2ANTI 12096
SMALTALG 12098
SME2SNSE 12098
SMMFLGTH 388
SME1ANTI 12097
SMSPLTRD 12099
SMBSMMRT 12098
SME1SNSE 12097
SME1PCTS 11725
SMRRNART 12097
SME1MPRT 11622
SMNUM5CD 1
SMDPMPRT 2
SME2PCTS 11604
```
### Some cleaning--- getting rid of all rows with nan values
```
df1 = df.iloc[:147,:]
df2 = df.iloc[160:7294,:]
df3 = df.iloc[8828:12365,:]
df4 = df.iloc[13375:13477,:]
df5 = df.iloc[13570:14674,:]
df6 = df.iloc[14768:14838,:]
df7 = df.iloc[14892:,:]
new_df = pd.concat([df1,df2,df3,df4,df5,df6,df7], axis = 0)
new_df.shape
(12102, 64)
```
### list od feautres to be droped from new_df dataframe
``` columns_to_drop = ['SMMNCPB','SMNUM5CD', 'SMGAPPCT', 'SMCGLGTH', 'SMMNCV','SMCENTER','SMTSPAX','SMGTC','SMNUMGPS','SM550NRM', 'SM350NRM', 'SMPTHNTS']
new_df.drop(columns = columns_to_drop, inplace = True)
new_df.shape
(12102, 52)
```
### It is always useful to get an overall idea of the whole dataset and narrow it down from there
### plotting the histograms of the whole numeric data
```
new_df.hist(figsize= (30,30), bins = 50)
plt.tight_layout
```
![Image_1](/img/13.png)
Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc



>>>>>>> 301a406916436cd9751b68343edc29a78123cff1
