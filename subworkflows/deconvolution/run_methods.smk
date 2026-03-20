
import os
import sys
import yaml

methods = config["methods"].split(',')   
# Inclure dynamiquement les pipelines appropriés en fonction des méthodes
for method in methods:
    method = method.strip()  
    Snakefile = f"{method}/run_{method}.smk"
    include: Snakefile