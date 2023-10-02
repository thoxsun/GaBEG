# GaBEG
Gaia bolometric Brightness Evaluation for Globular clusters 

For evaluation of the bolometric GC TRGB brightness with Gaia DR3 and future releases, I
wrote the Gaia bolometric Brightness Evaluation for Globular clusters (GaBEG) modules.
The GaBEG modules are made available under the GNU General Public License v3.0. 

The modules can be applied to any GC, and even to a
certain extend to OCs and other gravitational stellar bound systems. Note however, that the
modules are tested and optimized for the GC M5, and to a certain extend to GC M3 and OC
M67. Adjustments are necessary for the use on other systems. 


The GaBEG modules are dependent on: 
-astropy version V4.3.1 (Collaboration et al., 2022a), 
-astroquery V0.4.6 (Ginsburg et al., 2019), 
-scipy V1.5.3 (Virtanen et al., 2020), and
-matplotlib V3.5.2 (Caswell et al., 2022).



The GaBEG modules are separated in 
• cluster member selection modules, and 
• bolometric TRGB brightness evaluation modules, 
which depend on the group of fundamental modules.


**Cluster member selection modules**
Following cluster member selection modules are called the **single parameter selection
modules**, which can be used in arbitrary order:

• proper motion    (proper_motion_evaluation_loop.py)
• radial velocity  (gaiaRadV elSelection.py)
• parallax         (gaiaP arallaxSelection.py)

For **cluster member selection**, one most use the adjoined modules
1) angular position (GAIA_query_circle.py)
2) single parameter selection modules
3) cluster member selection (gaiaClusterM emberSelection.py)
in the listed order.


**Bolometric TRGB brightness evaluation modules**
To derive the bolometric TRGB brightness of a GC, the bolometric brightness evaluation modules

1) distance (gaiaDistEval.py)
2) extinction (gaiaExt.py)
3) bolometric correction with Gaiadr3 Bcg algorithm (gaiaBCg.py)
4) bolometric correction approximation (gaiaBCg_approximation.py)
5) RGB selection (gaiaAGBRGBseparation.py)
6) TRGB and axion-electron coupling bound evaluation (gaiaT RGB_g13_eval.py)
need to be executed in the listed order. 

The single object epoch photometry analysis module
• gaiaSingleObjectAnalysis.py 
enables the derivation of stars being variables and should be used to investigate the brightest RGB star, selected for the TRGB evaluation. 


The modules directly include documentation on necessary input data to load, package
dependencies, functionality, structure, parameters, and generated tables and plots.

Note: 
• The GaBEG modules are a prototype written for the evaluation of the results in my Master Thesis (Link will follow). For a detailed explanation on evaluation and the description of methods applied to the data, have a look at the thesis. All main graphs can be found there as well. 
• The modules might be updated in the future. For now, they remain a prototype and still need adjustment for generalization to any cluster. Moreover, documentation of some modules is sparse.
• Nonetheless, I hope that this modules (including KDE-based cluster member selection, interpolations, beautiful plots and more) will serve well as a basis and source for future evaluation of Gaia data. 
