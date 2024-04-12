# Rare-Event-Detection

This unsupervised machine learning algorithm learns to decompose a trajectory into a set of components that represent structural motifs that move simultaneously and reports their appearance as dominant structural characteristics during a particular period in the course of the trajectory. 

For trajectory analyses the input array (I) is a contact matrix (size: number_of_frames × number_of_pairs) constructed in step1-5 from the data in each frame of the trajectory. By construction, this array encodes the information about the protein structure as it evolves over the trajectory time. In step6, the NMF algorithm decomposes the I array into a “spatial” array (H) of size (number_of_components × number_of_pairs) and a “temporal” array (W) of size (number_of_frames × number_of_components), which together are responsible for the integrated detection of the temporally defined rare events that provide the mechanistic information.

The number_of_components is a user-defined parameter. The selection of the number_of_components should be guided by the number of independent processes in the system. A adequate number of components is evident when the processes of interest are effectively captured by the transitions of existing components, and adding more components would mainly capture fluctuations of lesser significance.

Moreover, I noticed an interplay between the choice of number_of_components and the choice of the regularization type parameter l1_ratio of the NMF step. With l1_ratio=1, the fitting is regularized by a L1 norm term, while with l1_ratio=0 the regularization is based on the Frobenius norm. Compared to L1 norm, the regularization based on the Frobenius norm tends to give a smaller penalty to more evenly distributed parameters. Thus, in RED analysis with l1_ratio=1, the results I obtained tends to identify a component that contributes consistently throughout the trajectory, which look like the result in ref(1). In contrast, with l1_ratio=0, this constitutive component do not emerge as a separate component. Instead, the constitutive structural elements are distributed across the other components. Thus, all components resulting from the RED analysis with l1_ratio=0 contain the same constitutive residue pairs, as in ref(2). Interestingly, while the result in ref(2) was obtained with number_of_components=5 and l1_ratio=0, if the NMF was applied on the same dataset with l1_ratio=1 that favors the appearence of a constitutive component and ith number_of_components=6, the result I got was a constitutive component plus 5 other components whose temperol evolution looks exactly the same as the result with number_of_components=5 and l1_ratio=0, and whose structural elements are also very similar except not having the constitutive pairs in every component. These comparison shows that the two types of RED analysis are equivalent in terms of the events they captured, with the difference lying only in the representation of the constitutive residue pairs. In the following analysis, I recommend using the version with l1_ratio=0 where the constitutive residue pairs are consistently appearing in every component, as it facilitates the data comparison and interpretation as described in ref(2).

References:
1. Plante, A., and Weinstein, H. (2021) Ligand-Dependent Conformational Transitions in Molecular Dynamics Trajectories of GPCRs Revealed by a New Machine Learning Rare Event Detection Protocol. Molecules 26
2. Xie, H., and Weinstein, H. (2023) Allosterically coupled conformational dynamics in solution prepare the sterol transfer protein StarD4 to release its cargo upon interaction with target membranes. Frontiers in Molecular Biosciences 10
