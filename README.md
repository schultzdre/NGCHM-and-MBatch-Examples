# NGCHM and MBatch Examples
This repository provides example code using the Next-Generation Clustered HeatMap ([NGCHM](http://www.ngchm.net/)) and MBatch R packages. To install the NGCHM R package follow the instructions at:

[https://github.com/bmbroom/NGCHMR](https://github.com/bmbroom/NGCHMR)

To install the MBatch package follow the instructions at:

[https://github.com/MD-Anderson-Bioinformatics/MBatch](https://github.com/MD-Anderson-Bioinformatics/MBatch)

## Data Generation
The code uses simulated data generated using the ```simulateData.R``` file. The data is simulated in two groups (e.g. cell types) and two batches (e.g. sequencing runs). The data is generated as follows:
1. Mean expression for each gene in group one are first generated from a Normal distribution.
2. Mean expression for each gene in group two are calculated by first generating random numbers from a normal distribution with mean one, then multiplying that number by the mean of the same gene in group one.
3. Gene expression values for each sample are then generated from a Normal distribution using the generated mean expression values.
4. Next, a number of samples are taken from each group and added to each batch to simulate technical replicates. The remaining samples are then distributed between batches according to the batch imbalance. Please refer to the examples section for an illustration of how samples are distributed.
5. Each numerical value is then multiplied by a random number to introduce noise.
6. The expression of each gene in each batch is then shifted by a random number sampled Normally and scaled by a random number sampled from and inverse gamma distribution to generate the final data.

The ```simulateData``` function has the following parameters:

**sampleNum** - Numeric. Number of samples generated.  
**replicateNum**- Numeric. Number of technical replicates taken per group. Final number of samples exported will be *sampleNum + 2\*replicateNum*.  
**batchImbalance** - Numeric. Number of samples added to one of the batches after accounting for technical replicates. All remaining samples from the group are added to the other batch.  
**geneNum** - Numeric. Number of genes to be simulated.  
**setTheSeed** - Numeric. Number used to set the random number generator seed. Default NULL.  
**groupOneMeanSD** - Numeric vector of length two. Mean and Standard Deviation respectively used to generate mean gene expression in group one. Default c(0,1.5).  
**groupTwoSD** - Numeric. Standard Deviation used to generate mean gene expression in group two. Default 1.  
**sampledValueSD** - Numeric. Standard deviation used to generate gene expression values using the calculated mean expression for each gene in each group. Default 1.  
**noiseLevel** - Numeric. Standard deviation of noise level. Sampled expression levels are multiplied by normal numbers taken from a distribution with mean one and standard deviation defined here. Default 0.25.  
**batchEffectLocation** - Numeric vector of length two. Mean and standard deviation respectively of location batch effects. Default c(1,2).  
**batchEffectScale** - Numeric Value of length four. Shape and scale of inverse gamma distributions used to sample scale batch effects. First two values identify shape and scale of distribution on batch one, and third and fourth values of batch two. Default c(6,4,4,5).

As an example, using the parameters *sampleNum* = 200, *replicateNum* = 15 and *batchImbalance* = 20, one hundred samples will be initially generated for each group. Next, 15 samples from each group will be used as technical replicates. These samples will be duplicated and added to each batch so each batch then has 30 replicate samples. Of the remaining 85 samples in group1, 20 are then added to batch1 and 65 to batch2. Similarly, of the remaining 85 samples in group2, 20 are added to batch2 and 65 to batch1.

## Running Example
The script ```runExamples.R``` generates a simulated dataset using ```sampleNum = 500```, ```replicateNum = 20```, ```batchImbalance = 100```, and ```geneNum = 400```. The example then corrects the data using both EBNplus and RBN. NGCHMs are plotted for the uncorrected and both corrected datasets. The ```.ngchm``` files can be viewed using the NGCHM viewer at [https://www.ngchm.net/Downloads/ngChmApp.html](https://www.ngchm.net/Downloads/ngChmApp.html).

The uncorrected datasets clusters first by batch then by group:
![](https://bioinformatics.mdanderson.org/Software/SAMMI/Thumbnails/uncorrected.PNG)

Using EBN correction we see that samples cluster by group, and technical replicates cluster next to each other:
![](https://bioinformatics.mdanderson.org/Software/SAMMI/Thumbnails/EBN_corrected.PNG)

Similar results are found for RBN correction:
![](https://bioinformatics.mdanderson.org/Software/SAMMI/Thumbnails/RBN_corrected.PNG)
