# sise.micro-rna-m-rna-targets-interactions.backend

*sise.micro-rna-m-rna-targets-interactions.backend*<br>
✅ [Visit our website here](http://132.73.84.177/).

Created By:<br>
 Uri Zlotkin<br>
 Shahar Kramer<br>

## Sections

- [Web site for miRNA-mRNA interactions data](#Web_site_for_miRNA-mRNA_interactions_data)
- [Getting Started](#Getting_Started)
- [Endpoints](#Endpoints)
- [Data Access Layer](#Data_Access_Layer)
- [Configuration](#Configuration)

## Web_site_for_miRNA-mRNA_interactions_data

The backend is developed with Python Flask, and the data is stored and handled with Postgres SQL.

## Getting_Started

In order to set up the development environment and run the project, we create a file called requirements.txt.
All the packages we use in the project are in that file.<br>
flask==2.2.2<br>
flask_cors==3.0.10<br>
pandas==1.5.2<br>
numpy==1.23.5<br>
flask-sqlalchemy<br>
psycopg2<br>
flask_executor<br>
flask-compress<br>
waitress

## Endpoints

132.73.84.177/api - welcome message for the api

132.73.84.177/api/organisms/details - return list of organisms in json.

132.73.84.177/api/organisms/datasets/<int:data_set_id>/interactions - return list of interactions of dataset:<int:data_set_id> in json

132.73.84.177/api/interactions - return list of interactions base on search filters:
(datasetsIds, seedFamilies, miRnaIds, miRnaSeqs, sites, geneIds, regions)           

## Data_Access_Layer

PostgreSQL is used as a database, and flask_sqlalchemy is used to communicate with it.
Each table is translated into a Python object by mapping it to a python object.
Our objects: (Organism, DataSet, Interaction, mirnaIdOption, SeedFamilyOption, GeneIdOption, RegionOption, SiteOption).

## Configuration

There is a config.ini file and a Configurator class.
Database connection information is stored in the config.ini file, which is not uploaded to git.
The Configurator class reads the necessary data from the config.in file.


