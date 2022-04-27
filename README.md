# GPU Acceleration of Set Similarity Joins

This repository contains the implementation of a scheme of efficient set similarity joins on Graphics Processing Units (GPUs). This solution takes advantage of the massive parallel processing offered by GPUs and employs MinHash to estimate the similarity between two sets in terms of Jaccard similarity. Experimental results show that our proposed method is more than two orders of magnitude faster than the serial version of CPU implementation, and 25 times faster than the parallel version of CPU implementation, while generating highly precise query results.

Related publications:

CRUZ, M. S. H.; AMAGASA, T.; WATANABE, C.; LU, W.; KITAGAWA, H.. Secure similarity joins using fully homomorphic encryption. In: Proceedings of the 19th International Conference on Information Integration and Web-based Applications & Services (IIWAS), Vienna – Austria, 2017 https://doi.org/10.1145/3151759.3151788

CRUZ, M. S. H.; KOZAWA, Y.; AMAGASA, T.; KITAGAWA, H.. Accelerating Set Similarity Joins Using GPUs. In: Transactions on Large-Scale Data and Knowledge-Centered Systems (TLDKS), 2016. http://link.springer.com/chapter/10.1007%2F978-3-662-53455-7_1

CRUZ, M. S. H.; KOZAWA, Y.; AMAGASA, T.; KITAGAWA, H.. GPU Acceleration of Set Similarity Joins. In: The 26th International Conference on Database and Expert Systems Applications (DEXA 2015), Valencia – Spain, 2015 http://link.springer.com/chapter/10.1007/978-3-319-22849-5_26
