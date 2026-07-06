
import os
import sys
import yaml

methods = config["methods"].split(',')   
# Dynamically include the appropriate pipelines based on the selected methods
for method in methods:
    method = method.strip()  
    Snakefile = f"{method}/run_{method}.smk"
    include: Snakefile