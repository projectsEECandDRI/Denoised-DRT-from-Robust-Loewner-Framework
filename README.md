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
