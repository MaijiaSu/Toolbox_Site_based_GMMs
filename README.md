# Toolbox: Site-based stochastic ground motion model (GMM) 
This project provides a toolbox to fit site-based stochastic ground motion models (GMMs) and generate statistically compatible synthetic ground motions (GMs). Site-based stochastic GMMs are fitted to a dataset of recorded GMs representing the seismic hazard of interest. The dataset can be selected with the following methods:
* **Traditional GM Selection Methods**: Records are selected to match target spectra, such as a uniform hazard spectrum or a conditional mean spectrum. The fitted stochastic GMM can then generate an unlimited number of synthetic records compatible with the chosen dataset.
* **Selection Based on Similar Seismological Parameters**: Records are selected based on specific earthquake characteristics (e.g., magnitude $M_w$, source-to-site distance $R$, and site condition $V_{S30}$), leveraging the ergodic assumption. The fitted GMM produces synthetic GMs that are statistically consistent with the selected records.
* **Selection from Ranges of Seismological Parameters**: Records are chosen within specified ranges of $M_w$, $R$, and $V_{S30}$. This method also relies on the ergodic assumption but requires developing vectorized GM prediction equations (GMPEs) to link GMM parameters to seismological variables, allowing scenario-specific GM simulations.

#  Brief introduction to site-based stochastic GMMs

#  How to use this toolbox?
## I. Data Collection
All selected ground motions are saved in a MATLAB structure variable, 
GM(ir).eq = 
GM(ir).dt =

MetaData


## II. Fit the GM Dataset
### Step 1： 
### Step 2: 
 * either Fiting a joint probablisitic density function
 * or Fiting a vectrozied GMPEs

## III. Generate synthetic GMs


## IV. Visualize the Results


# Cite our works



~~example~~,
$\mathcal{M}(\boldsymbol(\theta))$

- [ ] Task not done
- [x] Task completed
