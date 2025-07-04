# Toolbox: Site-based stochastic ground motion model (GMM) 
This project provides a toolbox for fitting site-based stochastic ground motion models (GMMs) and generating statistically compatible synthetic ground motions (GMs). Site-based stochastic GMMs are fitted to a dataset of recorded GMs representing the seismic hazard of interest. The dataset can be selected with the following methods:
* **Traditional GM Selection Methods**: Records are selected to match target spectra, such as a uniform hazard spectrum or a conditional mean spectrum. The fitted stochastic GMM can then generate an unlimited number of synthetic records compatible with the chosen dataset.
* **Selection Based on Similar Seismological Parameters**: Records are selected based on specific earthquake characteristics (e.g., magnitude $M_w$, source-to-site distance $R$, and site condition $V_{S30}$), leveraging the ergodic assumption. The fitted GMM produces synthetic GMs that are statistically consistent with the selected records.
* **Selection from Ranges of Seismological Parameters**: Records are chosen within specified ranges of $M_w$, $R$, and $V_{S30}$. This method also relies on the ergodic assumption but requires developing vectorized GM prediction equations (GMPEs) to link GMM parameters to seismological variables, allowing scenario-specific GM simulations.

#  How to use this toolbox?
## I. Data Collection
All selected ground motions are saved in two MATLAB structure variables, 

```matlab
for i = 1:NofGM % Number of Ground Motions
    for n = 1:NofComp % Number of Components
        GM{n,i}.eq = ... % acceleration signal
        GM{n,i}.dt = ... % time step
    end
    MetaData.M(i) = ... % Momment magnitude
    MetaData.VS30(i) = ... % Shear wave velocity
    MetaData.R(i) = ...    % source-to-site distance
    ...
end
```

## II. Fit the Selected Dataset
See  ```main_Fit.m ``` 

## III. Generate desirable synthetic ground motions
See  ```main_gen.m ``` 
