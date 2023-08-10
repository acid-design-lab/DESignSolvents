# DESignSolvents: Quantitative Prediction of Physicochemical Properties of Deep Eutectic Solvents Using Machine Learning

Deep eutectic solvents (DESs) have emerged as a promising eco-friendly alternative to environmentally hazardous organic solvents. DESs are a new class of mixtures characterized by a reduced melting point relative to pure starting materials. This phenomenon is explained by the formation of hydrogen bonds. However, selecting the suitable DES for a specific application has been traditionally challenging and time-consuming. This study presents an innovative approach to optimizing DES design by developing a machine learning workflow that enables us to predict fundamental physicochemical properties — melting point, density, and viscosity — critical to determining the scope of DES. Our models, particularly gradient-boosted trees like CatBoost, exhibited strong predictive accuracy, with cross-validation R2 values of 0.76, 0.89, and 0.64 for the mentioned properties, respectively. ore importantly, we created a web resource, DESignSolvents (ССЫЛКА), which hosts an extensive database of two- and three-component DESs properties and associated prediction models. This resource can significantly contribute to accelerating the development of DES for various applications, adopting green chemistry practices, and facilitating the global transition to sustainable solvent usage.

![image](https://github.com/Odegova-Valerie/DESignSolvents/assets/101416592/8d040edc-2d07-4f94-9d12-edcf28ff2480)


## Sections
1. DataBases
  Each folder contains databases containing the values of the physico-chemical properties of deep eutectic solvents. As well as databases for machine learning,   which contain all the descriptors used to predict the properties of DES

2. DataAnalysis
   In these notebooks, the primary processing of the database takes place: SMILES parsing, deletion of repetitions, as well as exploratory data analysis

3. MiningDescriptors
  This file contains code for getting descriptors using RDKit, as well as using various formulas

4. ML_MODELS
   This section contains code for optimizing machine learning models and further selecting the best model. This happens by counting various metrics such as R2 and RMSE

5. BEST_MODEL
  In these sections, the best model is presented, as well as checking it on test data. In addition, each folder contains the best trained model in pkl format


