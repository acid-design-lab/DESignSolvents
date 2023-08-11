# Quantitative Prediction of Physicochemical Properties of Deep Eutectic Solvents Using Machine Learning

Deep eutectic solvents (DESs) have emerged as a promising eco-friendly alternative to environmentally hazardous organic solvents. DESs are a new class of mixtures characterized by a reduced melting point relative to pure starting materials. This phenomenon is explained by the formation of hydrogen bonds. However, selecting the suitable DES for a specific application has been traditionally challenging and time-consuming. This study presents an innovative approach to optimizing DES design by developing a machine learning workflow that enables us to predict fundamental physicochemical properties — melting point, density, and viscosity — critical to determining the scope of DES. Our models, particularly gradient-boosted trees like CatBoost, exhibited strong predictive accuracy, with cross-validation R2 values of 0.76, 0.89, and 0.64 for the mentioned properties, respectively. ore importantly, we created a web resource, **DESignSolvents**, which hosts an extensive database of two- and three-component DESs properties and associated prediction models. This resource can significantly contribute to accelerating the development of DES for various applications, adopting green chemistry practices, and facilitating the global transition to sustainable solvent usage.


![image_2023-08-12_01-14-51](https://github.com/acid-design-lab/DESignSolvents/assets/82499756/e3df5679-abe7-4839-997a-47c3b48550f8)


## Repository contents

Each folder in this repository contains files specific to one of the three physico-chemical properties predicted in this work.

### Databases
`CSV` files are the database of the corresponding physico-chemical property and the dataset used to train machine learning models.

### Data engineering and analysis
`DataAnalysis` notebook contains the main processing steps, such as SMILES parsing, removal of duplicates, exploratory data analysis, etc.

### Mining descriptors
`MiningDescriptors` file contains the code for calculating descriptors using RDKit.

### ML models
`ML` notebook contains the code for optimizing and selectin the best machine learning model.

### Best model
`BEST_MODEL` notebook contains the optimized ML model of the top performance. The models are also available in the `PKL` format.


