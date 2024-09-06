
import os
import sys

methods = config["methods"].split(',')   
# Inclure dynamiquement les pipelines appropriés en fonction des méthodes
for method in methods:
    method = method.strip()  
    include: f"{method}/run_{method}.smk"