# HPLC-parameters-derivation-and-retention-time-prediction
Chromatography Retention Modeling: Quadratic and Neue–Kuss Implementation
This repository contains the source code and implementation details for predicting solute retention in Reverse-Phase High-Performance Liquid Chromatography (RP-HPLC). The project focuses on deriving and optimizing intrinsic retention parameters using two complementary modeling frameworks.
🔬 Scientific Background
The scripts implemented here facilitate the transition from isocratic experimental data to gradient retention prediction by solving the fundamental equation of gradient elution.
Supported Models
The Quadratic model is expressed as:

$$
\ln k(\varphi) = \ln k_{\text{w}} - S_1\varphi + S_2\varphi^2
$$

The Neue–Kuss model is expressed as:

$$
\ln k(\varphi) = \ln k_{\text{w}} + 2\ln(1 + c_{\text{nk}}\varphi) - \frac{b_{\text{nk}}\varphi}{1 + c_{\text{nk}}\varphi}
$$
                
💻Technical Implementation
All numerical optimizations and transport simulations were implemented using R (v4.4.2).            
Core MethodologyGradient Integration: To predict retention behavior under gradient conditions, the program solves the governing Ordinary Differential Equations (ODEs) describing solute transport.

Numerical Solver: The deSolve package is employed to perform robust integration of the solute migration velocity:

$$\frac{dz}{dt} = \frac{u}{1 + k(\varphi(t))}$$

where $u$ represents the mobile phase linear velocity and $k(\varphi(t))$ is the instantaneous retention factor. 

To account for the system dwell time ($t_d$), the effective mobile phase composition $\varphi(t)$ at the column head is modeled as a piecewise function:

$$\varphi(t) = 
\begin{cases} 
\varphi_0, & t \le t_d \\
\varphi_0 + B(t - t_d), & t > t_d 
\end{cases}$$

where $\varphi_0$ is the initial organic volume fraction and $B$ is the gradient slope.

Elution Condition: The elution point is reached when the cumulative migration distance equals the column length ($L$). The predicted retention time ($t_{\text{R}}$) is determined by solving the integral form of the migration equation:

$$\int_{0}^{t_{\text{R}}} \frac{u}{1 + k(\varphi(t))} \ dt = L$$

Optimization Strategy: Unlike conventional methods that minimize temporal residuals, parameter estimation in this study is performed by minimizing the Sum of Squared Errors (SSE) relative to the column length ($L$). This objective function ensures that the fundamental transport equation is satisfied at each experimentally observed retention time ($t_{R,\text{exp},i}$):

$$SSE = \sum_{i=1}^{n} \left( \int_{0}^{t_{R,\text{exp},i}} \frac{u}{1 + k(\varphi(t))} \ dt - L \right)^2$$

🛠️ Prerequisites
To run these scripts, you need R v4.4.2 or higher with the following packages installed:
install.packages(c("numDeriv","dplyr","parallel","readxl","writexl","deSolve"))

### 📂 File Structure
*   **[`/scripts`](./scripts)**: Contains the main R scripts for model fitting.
    *   [`Quadratic_Model_Fitting.R`](./scripts/Quadratic_Model_Fitting.R): Implementation of the Quadratic model.
    *   [`Neue_Kuss_Model_Fitting.R`](./scripts/Neue_Kuss_Model_Fitting.R): Implementation of the Neue–Kuss model.
*   **[`/data`](./data)**:
    *   [`example_data.xlsx`](./data/example_data.xlsx): Template dataset for chromatography runs.
*   **[`/results`](./results)**:
    *   `results.xlsx`: Auto-generated output file.

🚀 Usage
Data Preparation: Ensure your chromatographic data (retention times, gradient profiles, and column parameters) is formatted according to the template in data/example_data.xlsx.
Model Selection: Choose the appropriate script based on your research objective:
Use Quadratic_Model_Fitting.R for the second-order polynomial model.
Use Neue_Kuss_Model_Fitting.R for the Neue–Kuss model.
Execution: Open the selected script in RStudio or run it via the terminal:
Bash
Rscript scripts/Quadratic_Model_Fitting.R
Result Retrieval: Once the optimization converges, the calculated parameters ($\ln k_{\text{w}}$, $S$, etc.) will be exported to results/results.xlsx.
Note: The scripts utilize the parallel package to accelerate the iterative optimization process.


📝 Citation
If you use these scripts in your research, please cite this repository and the associated manuscript:

1. **Codebase:** Guan, L. (2026). *HPLC Retention Modeling Scripts*, v1.0.0, Zenodo, https://doi.org/10.5281/zenodo.xxxxxx
2. **Manuscript:** Guan, L., Han, C., Zhu, L., & Wang, B. (2026). *An R-Based Computational Framework Integrating Mechanistic LSER Modeling for Enhanced $\log P_{\text{ow}}$ Determination via Gradient Reversed-Phase HPLC*. (In preparation/Submitted).

📧 Contact
For technical inquiries regarding the R implementation or the ODE solver logic, please contact the corresponding author or the repository maintainer.


          
