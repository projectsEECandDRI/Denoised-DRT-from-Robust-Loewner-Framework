# Denoised-DRT-from-Robust-Loewner-Framework

This repository provides MATLAB scripts and a graphical user interface (GUI) for extracting denoised distribution of relaxation times (DRT) from noisy electrochemical impedance spectroscopy (EIS) data using a data-driven approach called the Robust Loewner Framework (RLF).

The RLF algorithm builds upon the direct interpolation capabilities of the Loewner Framework (LF) while addressing its sensitivity to noise by incorporating:
* **ScreeNOT**: A robust method for selecting the optimal SVD truncation threshold that determines the optimal model order for the reduced Loewner model.
* **Akaike Information Criterion (AIC)-based filtering**: A post-processing step that refines the optimal model order, yielding a physically meaningful, so-called denoised DRT.

The RLF method is well-suited for frequency response analysis and system characterization, both within electrochemical impedance spectroscopy (EIS) and beyond.

## 🔑 Key Features

- 🚫 No regularization or iterative optimization required 
- ⚡ Computationally efficient and accurate  
- 🙌 User-independent and parameter-free 
- 🔧 Versatile and easy to apply  

## 📁 Repository Contents

- Sample scripts to run the RLF algorithm on:
  - Noisy synthetic data from various Equivalent Circuit Models (ECMs)
  - Experimental data presented in the paper
- A user-friendly MATLAB GUI for interactive DRT extraction
- User manual with instructions for using the GUI
- Example datasets from:
  - Ferrocyanide oxidation experiments
  - PEM water electrolysis experiments

## 🛠️ Requirements

- MATLAB R2019b or later

## 📝 Algorithmic Steps for RLF-Based DRT Extraction from Impedance Data

### Step 1: Import Impedance Data

- Load the frequency vector $f$ and complex impedance data $Z(f) = Z\_{real} + i Z\_{imag}$

### Step 2: Construct Loewner Matrices

- Use `Loewner_Framework.m` to compute:
  - Loewner matrix $L$, shifted Loewner matrix $Ls$
  - Interpolation vectors $V$, $W$

### Step 3: Estimate Optimal Model Order for SVD Truncation

- Use `ScreeNOT.m` to determine optimal SVD truncation order $r\_{opt}$
  - Eliminates the need for manual thresholding
  - Enhances robustness to noise

### Step 4: Construct Reduced-Order State-Space Model

- Use `state_space_mod.m` to compute matrices ($E_{r\_{opt}}$, $A_{r\_{opt}}$, $B_{r\_{opt}}$, and $C_{r\_{opt}}$)
- These define the rational approximation:
  $Z_{r\_{opt}}(i\omega) = C_{r\_{opt}} (i\omega E_{r\_{opt}} - A_{r\_{opt}})^{-1} B_{r\_{opt}}$


### Step 5: Compute Preliminary DRT (Unfiltered)

- Use `DRT_values.m` to extract:
  - Relaxation times:  
    \[
    \tau_i = \frac{1}{\mathrm{Re}(\lambda_i)}
    \]  
    where \(\lambda_i\) are eigenvalues of \(E_r^{-1} A_r\)
  - Resistances \(R_i\) from eigen-decomposition
  - DRT set:  
    \[
    \{(\tau_i, R_i)\}_{i=1}^r
    \]

### Step 6: AIC-Based Model Filtering

- Use `AIC_calculation.m` to refine the DRT:
  - For each model order \(k \leq r\), compute:  
    \[
    \mathrm{AIC}(k) = n \cdot \log\left(\frac{\mathrm{RSS}_k}{n}\right) + 2k
    \]  
    where RSS is residual sum of squares, and \(n\) is number of frequency points
  - Select model order \(k^*\) minimizing AIC
  - Retain only top \(k^*\) DRT elements → denoised DRT

### Step 7: Reconstruct Impedance from Denoised DRT

- Approximate impedance using:  
  \[
  Z_{\text{DRT}}(s) = R_{\infty} + \sum_{i=1}^{k^*} \frac{R_i}{1 + s \tau_i}
  \]
- (Optional) Estimate \(R_{\infty}\) from high-frequency asymptote

### Step 8: Compute Relative Residual

- Evaluate reconstruction accuracy:  
  \[
  \text{Relative Residual} = \frac{\| Z_{\text{exp}} - Z_{\text{DRT}} \|_2}{\| Z_{\text{exp}} \|_2}
  \]
- Quantifies fidelity of extracted model to original data

### Step 9: Visualization and Analysis

- Plot:
  - Nyquist and Bode plots: compare original \(Z(f)\) and reconstructed \(Z_{\text{DRT}}(f)\)
  - DRT plot: \(|R_i|\) vs. \(\log_{10}(\tau_i)\)


## 📄 Citation

This repository is associated with the following research article:

**Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, Athanasios C. Antoulas, Tanja Vidaković-Koch**  
**A Data-driven, Noise-resilient Algorithm for Extraction of Distribution of Relaxation Times Using the Loewner Framework**  
**Journal of Power Sources, 2025.**

## 📚 Related Publications

If you use this repository or build upon this work, please also consider citing:

1. **Bansidhar Patel, Antonio Sorrentino, and Tanja Vidaković-Koch**,  
   *Data-driven analysis of electrochemical impedance spectroscopy using the Loewner framework*,  
   *iScience*, Vol. 28, No. 3, p. 111987, 2025.  
   [https://doi.org/10.1016/j.isci.2025.111987](https://doi.org/10.1016/j.isci.2025.111987)

2. **Antonio Sorrentino, Bansidhar Patel, Ion Victor Gosea, Athanasios C. Antoulas, Tanja Vidaković-Koch**,  
   *Determination of the distribution of relaxation times through Loewner framework: A direct and versatile approach*,  
   *Journal of Power Sources*, Vol. 585, 2023, 233575.  
   [https://doi.org/10.1016/j.jpowsour.2023.233575](https://doi.org/10.1016/j.jpowsour.2023.233575)

3. **Antonio Sorrentino, Ion Victor Gosea, Bansidhar Patel, Athanasios Antoulas, and Tanja Vidaković-Koch**,  
   *Loewner framework and distribution of relaxation times of electrochemical systems: Solving issues through a data-driven modeling approach*,  
   SSRN, 2022.  
   [https://ssrn.com/abstract=4217752](https://ssrn.com/abstract=4217752)  
   [http://dx.doi.org/10.2139/ssrn.4217752](http://dx.doi.org/10.2139/ssrn.4217752)
   
4. **Bansidhar Patel**,  
   *Application of Loewner framework for data-driven modeling and diagnosis of polymer electrolyte membrane fuel cells*,  
   *Master’s Thesis*, Otto-von-Guericke-Universität Magdeburg, 2021.  
   [https://doi.org/10.25673/101328](https://doi.org/10.25673/101328)

If you have questions, feedback, or need assistance, feel free to contact:  
📧 **patel@mpi-magdeburg.mpg.de**

In case you encounter any bugs or issues related to the GUI, please report them via email or open an issue on this repository.
