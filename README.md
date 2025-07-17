# Denoised-DRT-from-Robust-Loewner-Framework

This repository provides MATLAB scripts and a graphical user interface (GUI) ***(to be released soon)*** for extracting denoised distribution of relaxation times (DRT) from noisy electrochemical impedance spectroscopy (EIS) data using a data-driven approach called the Robust Loewner Framework (RLF).

This repository supports and accompanies the research work presented in the paper: 
**“A Data-driven, Noise-resilient Algorithm for Extraction of Distribution of Relaxation Times Using the Loewner Framework,”**  
*Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, Athanasios C. Antoulas, and Tanja Vidaković-Koch,*  
*Journal of Power Sources, vol. xxx, Art. no. 237909, 2025.*  
[https://doi.org/10.1016/j.jpowsour.2025.237909](https://doi.org/10.1016/j.jpowsour.2025.237909)

---

## 🚀Overview

The RLF algorithm builds upon the direct interpolation capabilities of the Loewner Framework (LF) while addressing its sensitivity to noise by incorporating:
* **ScreeNOT**: A robust method for selecting the optimal SVD truncation threshold that determines the optimal model order for the reduced Loewner model.
* **Akaike Information Criterion (AIC)-based filtering**: A post-processing step that refines the optimal model order, yielding a physically meaningful, so-called denoised DRT.

The RLF method is well-suited for frequency response analysis and system characterization, both within electrochemical impedance spectroscopy (EIS) and beyond.

---

## 🔑 Key Features

- 🚫 No regularization or iterative optimization required 
- ⚡ Computationally efficient and accurate  
- 🙌 User-independent and parameter-free 
- 🔧 Versatile and easy to apply  

---

## 📁 Repository Contents

- **Sample MATLAB scripts** for applying the RLF algorithm on:
  - Noisy synthetic data from various Equivalent Circuit Models (ECMs) discussed in the paper
  - Experimental data presented in the paper
- A user-friendly **MATLAB GUI** for interactive DRT extraction ***(to be released soon)***
- **User manual** with instructions for using the GUI ***(to be released soon)***
- **Example datasets** including:
  - Sample synthetic noisy EIS data from various ECMs discussed in the paper
  - Ferrocyanide oxidation experimental data from *(Patel, B. et al., 2025, [DOI](https://doi.org/10.17617/3.MMLZAK), Edmond V2)*  
  - PEM water electrolysis experimental data from *(Milicic, T. et al., 2024, [DOI](https://doi.org/10.17617/3.I0ECY9), Edmond V1)*

---

## 🛠️ Requirements

- MATLAB R2019b or later

---

## 📝 Algorithmic Steps for RLF-Based DRT Extraction from Impedance Data
The script **analyze_your_data_with_RLF.m** provides an efficient and user-friendly pipeline to extract a denoised DRT from your impedance data using the RLF approach. 
It performs the following main steps:
### Step 1: Import Impedance Data

- Load the frequency vector $f$ and complex impedance data $Z(f) = Z\_{real} + i Z\_{imag}$

### Step 2: Construct Loewner Matrices

- Use `Loewner_Framework.m` to compute:
  - Loewner matrix $L$, shifted Loewner matrix $Ls$
  - Interpolation vectors $V$, $W$

### Step 3: Estimate Optimal Model Order for SVD Truncation

- Use `ScreeNOT.m` (see [Acknowledgements](#acknowledgements) for source) to determine optimal SVD truncation order $r\_{opt}$
  - Eliminates the need for manual thresholding
  - Enhances robustness to noise

### Step 4: Construct Reduced-Order State-Space Model

- Use `state_space_mod.m` to compute matrices ($E_{r\_{opt}}$, $A_{r\_{opt}}$, $B_{r\_{opt}}$, and $C_{r\_{opt}}$)
- These define the rational approximation:
  $Z_{r\_{opt}}(i\omega) = C_{r\_{opt}} (i\omega E_{r\_{opt}} - A_{r\_{opt}})^{-1} B_{r\_{opt}}$

### Step 5: Extract Discrete DRT (Unfiltered)

- Use `DRT_values.m` to compute relaxation times ($\tau_i$) and corresponding resistances ($R_i$) by performing eigenvalue decomposition (EVD) of $-A_{r_{\text{opt}}}^{-1} E_{r_{\text{opt}}}$
- Plotting $|R_i|$ vs. $\log_{10}(\tau_i)$ yields the discrete DRT

### Step 6: AIC-Based Filtering

- First, sort the DRT components ($R_i$ and their corresponding $\tau_i$) in descending order based on the 2-norm of their residues. This step prioritizes poles ($-1/\tau_i$) that contribute most significantly to the system response, as indicated by larger residue norms:

$${||Res_i||}_2 = {||\frac{R_i}{\tau_i}||}_2$$

- Next, apply the corrected Akaike Information Criterion (AICc) using the script `AIC_calculation.m` to identify and retain only the most essential DRT components (referred to as the Denoised DRT).
- For each model order ( $r \le r_{\text{opt}}$), compute the corrected AIC as:

$${AICc(r) = N\cdot\log\left(\frac{{RSS}_r}{N}\right) + \frac{2r(r+1)}{N - r - 1} + 2r}$$
  
- where:   
  - $r$ is the number of parameters (i.e., the model order),  
  - $N$ is the number of observations (frequency points),  
  - $RSS_r$ is the residual sum of squares for model order ($r$ ).

- Determine the refined optimal model order ($r_{\text{opt}}^*$) by either:
  - selecting the model with the lowest AICc value, or  
  - identifying the largest drop in AICc between successive model orders.

- Retain the top $r_{\text{opt}}^*$ components to obtain the denoised DRT.


### Step 7: Reconstruct Impedance from Denoised DRT

- Approximate the impedance as:

$$Z_{\text{DRT}}(i\omega) = \sum_{i=1}^{r_{\text{opt}}^*} \frac{R_i}{1 + i\omega \tau_i}$$

- Alternatively, explicitly separate the ohmic resistance ($R_0$):
  
$$Z_{\text{DRT}}(i\omega) = R_0 + \sum_{i=1}^{r_{\text{opt}}^*-1} \frac{R_i}{1 + i\omega \tau_i}$$

### Step 8: Compute Relative Residual

- Assess the reconstruction accuracy via:

$$\text{Relative Residual} = \frac{{|| Z - Z_{DRT} ||}_2}{{|| Z ||}_2} $$

### Step 9: Visualization and Interpretation

- Generate key diagnostic plots for thorough analysis:
  - **Nyquist plot**: compare measured impedance **$Z(f)$** with reconstructed impedance **$Z_{\text{DRT}}(f)$**.
  - **DRT plot**: Visualize the denoised discrete DRT (**$|\mathbf{R_i}|$** vs. **$\log_{10}(\mathbf{\tau_i})$**).
  - **Relative residuals**: Quantify the discrepancy between measured and reconstructed impedance.
  - **Singular value decay**: Examine the decay of singular values from the Loewner pencil.
  - **AICc trends**: Visually inspect the corrected AIC values to confirm the selection of the refined optimal model order ($r_{\text{opt}}^*$).

---

## 📄 Citation

If you use this repository, please cite the following article:

**Bansidhar Patel, Antonio Sorrentino, Ion Victor Gosea, Athanasios C. Antoulas, Tanja Vidaković-Koch**  
*A Data-driven, Noise-resilient Algorithm for Extraction of Distribution of Relaxation Times Using the Loewner Framework,*  
*Journal of Power Sources, vol. xxx, Art. no. 237909, 2025.*  
[https://doi.org/10.1016/j.jpowsour.2025.237909](https://doi.org/10.1016/j.jpowsour.2025.237909)

---

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

---

## 🙏 Acknowledgements

The `ScreeNOT.m` script included in this repository is adapted from the MATLAB implementation available at [https://github.com/eladromanov/ScreeNOT](https://github.com/eladromanov/ScreeNOT). This implementation is based on the work by Donoho, D., Gavish, M., and Romanov, E. (2023), *"ScreeNOT: Exact MSE-optimal singular value thresholding in correlated noise,"* published in *The Annals of Statistics*, **51**(1), 122–148. We gratefully acknowledge the authors for making their method and code publicly available.

---

## 📬 Contact

For questions, feedback, or support, please contact:  
**Bansidhar Patel** — 📧 [patel@mpi-magdeburg.mpg.de](mailto:patel@mpi-magdeburg.mpg.de)

If you encounter bugs or issues with the GUI, please report them by email or open an issue in this repository.
