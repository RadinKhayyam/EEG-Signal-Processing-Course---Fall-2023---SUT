# 📊 EEG Signal Processing


**🏛️ University:** Sharif University of Technology  
**🏢 Department:** Electrical Engineering  
**👩‍🏫 Instructor:** Prof. Sepideh Hajipour  

---

## 📘 Course Overview  
This course provides an in-depth understanding of **EEG (Electroencephalography) signal processing**, emphasizing both the theoretical foundations and practical applications. Students will learn how to handle EEG data, preprocess signals, apply advanced source separation and classification techniques, and analyze brain connectivity. The course also covers the development of Brain-Computer Interfaces (BCI) and brain source localization methods. Practical projects and case studies are integrated into the curriculum to ensure hands-on experience in the field.

**Textbooks:**
- Sanei, Saeid. “Adaptive processing of brain signals.” John Wiley & Sons, 2013.
- Cohen, Mike X. Analyzing neural time series data: theory and practice. MIT press, 2014.
- Clerc, Maureen, Laurent Bougrain, and Fabien Lotte, eds. Brain-Computer Interfaces 1: Methods and Perspectives. John Wiley & Sons, 2016.
- Clerc, Maureen, Laurent Bougrain, and Fabien Lotte, eds. Brain-Computer Interfaces 2: Technology and Applications. John Wiley & Sons, 2016.

---

## 📝 Syllabus

### 1️⃣ Introduction to EEG Signal Processing  
- Overview of EEG signals, fundamentals of neuroelectric activity, and challenges in processing.  
- Introduction to key applications such as neurofeedback, BCI systems, and clinical diagnostics.

### 2️⃣ Functional Neuroimaging Techniques  
- Techniques for imaging brain activity, including comparisons of EEG with other modalities like fMRI and MEG.  
- The role of EEG in neuroimaging, its advantages, and limitations.

### 3️⃣ EEG Signals Pre-Processing  
- **Artifacts in EEG:** Overview of various artifacts (e.g., EOG, EMG, ECG) and their appearance in time-domain signals.  
- **Rereferencing Methods:** Common techniques for improving signal quality.  
- **Filtering:** Notch filter, bandpass filters, and others.  
- **Trial Rejection:** Methods for identifying and removing poor trials.  
- **Artifact Removal Using ICA:** Applying ICA to eliminate EEG artifacts by analyzing components across time-domain signals, frequency spectra, and spatial topography.

### 4️⃣ Blind and Semi-Blind Source Separation Methods  
- **Mathematical Basis of Blind Methods:** Principal Component Analysis (PCA) and Independent Component Analysis (ICA).  
- **Semi-Blind Methods:** Generalized Eigenvalue Decomposition (GEVD) and Denoising Source Separation (DSS).  
- Application of these methods for EEG noise reduction, comparing their effectiveness and evaluating noise reduction performance.

### 5️⃣ Processing of EEG Signals for Classification Purposes  
- **Feature Extraction:** Methods to extract relevant information from EEG data, such as power spectral density and time-domain features.  
- **Dimensionality Reduction:** Principal Component Analysis (PCA), Fisher’s Linear Discriminant (FLD), and various search-based methods.  
- **Classification Algorithms:** Simple classifiers like K-Nearest Neighbors (KNN), Support Vector Machines (SVM), and Naive Bayes.  
- **Evaluation Metrics:** Accuracy, precision, recall, F1-score, and ROC curves for assessing classification performance.

### 6️⃣ Patterns of Brain Signals and Detection Methods  
- **EEG Bands:** Delta, Theta, Alpha, Beta, and Gamma bands and their significance.  
- **Event-Related Potentials (ERP):** Introduction to the P300 wave and its detection methods.  
- **Steady-State Visual Evoked Potentials (SSVEP):** Classical detection methods and Canonical Correlation Analysis (CCA)-based approaches.  
- **Event-Related Synchronization (ERS):** Detection techniques for synchronization patterns.  
- **Feature Extraction Methods:** Using Common Spatial Patterns (CSP) for pattern recognition in EEG signals.

### 7️⃣ Brain-Computer Interfaces (BCI)  
- **Review of BCI Systems:** In-depth review of several articles on BCI, with a focus on P300-Spellers and SSVEP-Spellers.  
- **BCI Design and Implementation:** Steps in creating a BCI system, from signal acquisition to classification and feedback.

### 8️⃣ Brain Source Localization  
- **Forward and Inverse Problems:** Introduction to forward and inverse problems in brain source localization.  
- **Non-Parametric Solutions for Inverse Problem:** Minimum Norm Estimate (MNE), Weighted Minimum Norm Estimate (WMNE), Low-Resolution Electromagnetic Tomography (LORETA), and Standardized LORETA (sLORETA).  
- **Parametric Solutions for Inverse Problem:** Various parametric techniques for solving the inverse problem, including the investigation of the Mean Squared approach and introduction to Multiple Signal Classification (MUSIC)-based methods.  
- **Performance Comparison:** Comparing the performance and accuracy of different localization methods.

### 9️⃣ Brain Connectivity  
- **Types of Connectivity:** Structural, functional, and effective connectivity.  
- **Challenges in EEG Connectivity Analysis:** Addressing issues like volume conduction and active common reference effects.  
- **Functional Connectivity Estimators:** Pearson correlation, mutual information, coherence, Phase-Locking Value (PLV), Phase Lag Index (PLI), etc.  
- **Brain Connectivity Graphs:** Constructing connectivity graphs and extracting graph-theoretical features such as nodal degree, clustering coefficient, global efficiency, etc.

---
