# Supplementary data and R-scripts
Supplementary material for "Evaluating the predicted extinction risk of living amphibian species with the fossil record" [insert link to paper]

**This repository will be filled with the datafiles on paper acceptance.**

Includes datafiles, bib-file and Rmd-files to recreate the analysis and supplement file, and the supplement file itself as pdf (in case compilation fails).

## File descriptions
| File name                 | Description                        |
|---------------------------|------------------------------------|
| abbreviations.json        | Json file containing journal abbreviations to complement the csl file|
| bodysize.csv              | Body size data for extinct species |
| bodysize_caudata.csv      | Body size data for living caudate species collected from AmphibiaWeb|
| comb2.csv                 | Occurrence level data on species names and taxonomy, collection and occurrence numbers, lithology information, geological stage, stratigraphic range, (paleo-)coordinates and abundance counts |
| coords_from_shapefiles.RData | All coordinates extracted from shapefiles of geographic ranges of living amphibians (shapefiles retrieved from IUCN Red List webpage) |
| ecology_letters.csl       | Style file for *Ecology letter* citation style |
| ELEtietjeSA1.pdf          | Supplementary Information (figures, tables, modelling output, method part on data imputation |
| ELEtietjeSA1.Rmd          | R Markdown file to create the supplementary information pdf-file |
| European_Amphibians_Database3.csv | Body size data for living amphibian species. Supplement file from Trochet et al. 2014 [https://doi.org/10.3897/BDJ.2.e4123] |
| flowchart_model_building_and_cross_validation.jpg | Flowchart on model building and cross correlation procedure |
| Hirschfeld_Roedel_bodysize_data.csv | Body size data for living amphibian species. Supplement file from Hirschfeld & Roedel 2017: https://doi.org/10.1186/s12898-017-0135-y |
| imputed_data_extinct.RData | Imputation data from mice- and rF-algorithm approaches|
| iucn_abundance_scraping.RData | A presence abscence matrix for abundance keywords from text mining the population descriptions on the IUCN Red List webpage |
| iucn_export-amphibia-03feb2017.csv | IUCN Red List data on amphibian species including taxonomy, Red List status and assesssment details |
| living_data.csv           | Living species trait dataset that is used for prediction |
| model_data.csv            | Extinct species trait dataset that is used for model building and cross-validation |
| MS2.bib                   | bibliography file containing all references used in the manuscript and supplementary files |
| null-boot.RData           | Results from bootstrap on the Null model. As long as this file is present in the working directory the bootstrap is not performed again in Tietje_Roedel_2017_model_building_and_prediction.R as it lasts for about 20 minutes |
| occurrence_map.png        | Map of all fossil occurrences |
| Ruland-Jeschke_bodysize_data.csv |  Body size data for living amphibian species. Supplement file from Ruland & Jeschke 2016: https://doi.org/10.1016/j.biocon.2016.11.027 |
| stratplot.jpg             | Plots of stratigraphic range distributions |
| taxonomy.csv              | Dataframe with taxonomic information for extinct species |
| Tietje_Roedel_2017_data_processing.R | R script for data processing of extinct and living species traits |
| Tietje_Roedel_2017_model_building_and_prediction.R | R script for data imputation, model building, cross validation and prediction |
| Tietje_Roedel_2017_model_building_and_prediction_workspace.RData | R workspace created from the similarly named R script. This workspace is required to run ELEtietjeSA1.Rmd |
| Tietje_Roedel_2017_pbdb_references.csv | References from entries in the pbdb that were used in the final dataset |
| Tietje_Roedel_2017_trait_refs.pdf | All references for data collected literature or other databases |
