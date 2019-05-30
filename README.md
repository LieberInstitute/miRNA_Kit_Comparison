# miRNA_Kit_Comparison
This project evaluates 4 different commercially available small RNA sequencing kits and two computational preprocessing methods, one to collapse data to the reads present before amplification and a random subset control to determine if differences between the deduped data and the raw data was simply due to using fewer reads.


There are 4 main analyses:
1.	Similarity

    This was evaluated by using the same synthetic data and the same human brain data and comparing the miRNA estimates from the different methods.
    
2.	Accuracy

This was evaluated using equimolar synthetic oligonucleotides that are the same sequences as canonical human, rat, mouse, and virus miRNAs. We assessed how similar the results were across the different oligonucleotides.

3.	Detection Diversity

This was evaluated by testing the number of miRNAs and isomiRs detected by each method for both the synthetic miRNA and Brain samples.

4.	Consistency

    1.	Within a single batch and across triplicates
    2.	Across different batches
    
This was evaluated by looking at the correlation of counts for the same sample multiple times within the same batch or in different batches.


The preprocessing for these analysis was all done with shell scripts and the statistical analysis was performed using R.
See the **Final_Code** Directory for all code.
