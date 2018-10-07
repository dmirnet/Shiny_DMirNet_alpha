# Shiny_DMirNet_alpha
Shiny-based Web application for exploring direct microRNA-mRNA associations from expression profiles that supports parallel processing on cluster of computers
## Getting Started
These instructions will get you a copy of the web application up and running on your machines for development and testing purposes.
### Prerequisites
The following application must be installed to run DMirNet.
* [R Software](https://cran.r-project.org/) on each machines
* The following R packages must be install on the master machine to run Shiny_DMirNet. Run the following scripts to install the packages.
```R
install.packages("shiny")
install.packages("checkpoint")
```
## Running Shiny_DMirNet
There are two ways of installing DMirNet
* By running Script in R studio: To run DMirNet locally from GitHub, run the following script
```R
shiny::runGitHub('Shiny_DMirNet_alpha','dmirnet')
```
* By downloading the source code: Download the zip file of Shiny_DMirNet and extract the source file. Then run global.R, ui.R, or server. R source files. To run the source file, run the following script
```R
source("<path of the source file>/global.R")
```

## Citation
DMirNet is provided under a free-of-charge, open-source license (A-GPL3). All we require is that you cite/attribute the following in any work that benefits from this code or application.
### Citing the Web Application
Minsu Lee, Anene Bekuma Terefe, HyungJune Lee, and Sangyoon Oh, "DMirNet: Shiny-based application for exploring direct microRNA-mRNA associations", Bioinformatics (Under Review)
### Citing the DMirNet framework 
Minsu Lee and HyungJune Lee, (2016) "DMirNet: Inferring direct microRNA-mRNA association networks", BMC Systems Biology, 10(suppl5):125. 
