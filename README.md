# opm
This repository includes the code of an algorithm, IGV, proposed for the OpM problem, published in the following paper.

Mousavi, S., Bhambar, S., & England, M. (2023). An iterated greedy algorithm with variable reconstruction size for the obnoxious p‐median problem,  International Transactions in Operational Research, https://doi.org/10.1111/itor.13340

The algorithm is implemented in Java, C++, and Python. Only the last two languages were used in the paper though. Please adjust the algorithm based on your requirements (as we did differently in different sections of the paper). In particular, please adjust the termination condition. Currently, it is such that the algorithm runs 10 seconds on each instance, excluding the time needed to read the input data. Also, please change the path to where your data are stored and generate the report as you require. Currently, it only reports the best objective value and the duration (i.e., the running) time for each run on each instance. It currently runs 10 times per instance, which you may change. Finally, it currently runs on the large dataset described in the paper. 

In Section 5.1 of the paper, the original source of the small dataset and the original source of the pmed21.txt to pmed40.txt files, which were used to generate the large dataset, are mentioned. However, for convenience, these are also made available here, in the Data.zip file. 


Any comments and questions are welcome via email to seyed_r_mousavi@yahoo.com, as long as this email account remains active. 

Thank you.

