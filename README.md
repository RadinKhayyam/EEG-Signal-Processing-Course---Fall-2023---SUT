# EEG Signal Processing

**Department:** Electrical Engineering  
**University:** Sharif University of Technology  
**Instructor:** Prof. Sepideh Hajipour  

## Course Overview  
This course provides an in-depth understanding of EEG (Electroencephalography) signal processing, emphasizing both the theoretical foundations and practical applications. Students will learn how to handle EEG data, preprocess signals, apply advanced source separation and classification techniques, and analyze brain connectivity. The course also covers the development of Brain-Computer Interfaces (BCI) and brain source localization methods. Practical projects and case studies are integrated into the curriculum to ensure hands-on experience in the field.

## Syllabus

### 1. Introduction to EEG Signal Processing  
   - Overview of EEG signals, fundamentals of neuroelectric activity, and challenges in processing.  
   - Introduction to key applications such as neurofeedback, BCI systems, and clinical diagnostics.

### 2. Functional Neuroimaging Techniques  
   - Techniques for imaging brain activity, including comparisons of EEG with other modalities like fMRI and MEG.  
   - The role of EEG in neuroimaging, its advantages, and limitations.

### 3. EEG Signals Pre-Processing  
   - Filtering, artifact removal, and preparation of raw EEG data for analysis.  
   - Techniques such as notch filtering, bandpass filtering, and independent component analysis (ICA) for noise reduction.

### 4. Blind and Semi-Blind Source Separation Methods  
   - **Mathematical Basis of Blind Methods:** Principal Component Analysis (PCA) and Independent Component Analysis (ICA).  
   - **Semi-Blind Methods:** Generalized Eigenvalue Decomposition (GEVD) and Denoising Source Separation (DSS).  
   - Application of these methods for EEG noise reduction, comparing their effectiveness and evaluating noise reduction performance.

### 5. Processing of EEG Signals for Classification Purposes  
   - **Feature Extraction:** Methods to extract relevant information from EEG data, such as power spectral density and time-domain features.  
   - **Dimensionality Reduction:** Principal Component Analysis (PCA), Fisher’s Linear Discriminant (FLD), and various search-based methods.  
   - **Classification Algorithms:** Simple classifiers like K-Nearest Neighbors (KNN), Support Vector Machines (SVM), and Naive Bayes.  
   - **Evaluation Metrics:** Accuracy, precision, recall, F1-score, and ROC curves for assessing classification performance.

### 6. Patterns of Brain Signals and Detection Methods  
   - **EEG Bands:** Delta, Theta, Alpha, Beta, and Gamma bands and their significance.  
   - **Event-Related Potentials (ERP):** Detection of P300 components.  
   - **Steady-State Visual Evoked Potentials (SSVEP):** Classical detection methods and Canonical Correlation Analysis (CCA)-based approaches.  
   - **Event-Related Synchronization (ERS):** Detection techniques for synchronization patterns.  
   - **Feature Extraction Methods:** Using Common Spatial Patterns (CSP) for pattern recognition in EEG signals.

### 7. Brain-Computer Interfaces (BCI)  
   - **Review of BCI Systems:** In-depth review of several articles on BCI, with a focus on P300-Spellers and SSVEP-Spellers.  
   - **BCI Design and Implementation:** Steps in creating a BCI system, from signal acquisition to classification and feedback.

### 8. Brain Source Localization  
   - **Forward and Inverse Problems:** Mathematical formulations of forward and inverse problems in EEG.  
   - **Non-Parametric Methods:** Minimum Norm Estimate (MNE), Weighted Minimum Norm Estimate (WMNE), Low-Resolution Electromagnetic Tomography (LORETA), and Standardized LORETA (sLORETA).  
   - **Parametric Methods:** Various parametric techniques for solving the inverse problem.  
   - **Performance Comparison:** Comparing the performance and accuracy of different localization methods.

### 9. Brain Connectivity  
   - **Types of Connectivity:** Structural, functional, and effective connectivity.  
   - **Challenges in EEG Connectivity Analysis:** Addressing issues like volume conduction and active common reference effects.  
   - **Functional Connectivity Estimators:** Pearson correlation, mutual information, coherence, Phase-Locking Value (PLV), Phase Lag Index (PLI), etc.  
   - **Brain Connectivity Graphs:** Constructing connectivity graphs and extracting graph-theoretical features such as nodal degree, clustering coefficient, and global efficiency.

---

This syllabus reflects the graduate-level depth and rigor of the EEG Signal Processing course. The repository includes code and projects corresponding to each section of the syllabus.
