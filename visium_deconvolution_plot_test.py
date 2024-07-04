#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 16:46:25 2024

@author: maialen

Script to generate an interactive plot to show the top 5 cell types obtained
from an spatial deconvolution analysis

Input: 
    => CSV file containing normalized weights [barcodes as rownames, cell types as colnames]
    => CSV file containing spatial coordinates [barcodes as rownames, x and y coordinates as columns]

"""
# <! ------------------------------------------------------------------------!>
# <!                       IMPORTS                                           !>
# <! ------------------------------------------------------------------------!>
import pandas as pd
from math import pi
from bokeh.plotting import figure, output_file, save
from bokeh.transform import cumsum
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.palettes import Category20


# <! ------------------------------------------------------------------------!>
# <!                           DATA PREPARATION                              !>
# <! ------------------------------------------------------------------------!>

# Read spatial deconvolution result CSV file
norm_weights_filepath = "/home/maialen/Desktop/WORKSPACE/DECONVOLUTION/RCTD_DECONVOLUTION/GBM_UKF242_T_ST/norm_celltype_weights_UKF242_T_ST.csv"
norm_weights_df = pd.read_csv(norm_weights_filepath).set_index('Unnamed: 0')
norm_weights_df.index.name = None
#print(norm_weights_df.head())

# Read spatial coordinates CSV file
st_coords_filepath = "/home/maialen/Desktop/WORKSPACE/DECONVOLUTION/RCTD_DECONVOLUTION/GBM_UKF242_T_ST/spatial_coords_UKF242_T_ST.csv"
st_coords_df = pd.read_csv(st_coords_filepath).set_index('Unnamed: 0')
st_coords_df.index.name = None
#print(st_coords_df.head())

# Merge coordinate df and cell weight df
merged_df = pd.concat([st_coords_df, norm_weights_df], axis = 1)

# It will be difficult to show the information of all 54 cell types when hovering
# Thus, for each barcoded spot, retrieve the maximum 5 weights and create new columns
# accordingly. Those 5 max columns will be the info shown in the hovertool

# Get largest 5 cell types
n_largest_cell_types = 5
max_weights = norm_weights_df.apply(lambda x: x.nlargest(n_largest_cell_types).index.values, axis=1)


# Create df columns with max cell types
cell_type_storage_arrays = list()
cell_value_storage_arrays = list()

# Extract cell types with largest weights
for i in range(n_largest_cell_types):
    
    cell_type_storage_array = list()
    cell_value_storage_array = list()
    
    for barcode in max_weights.index:
        
        max_cell_types = max_weights.loc[barcode]
        max_cell_type = max_cell_types[i]
        max_cell_value = merged_df.loc[barcode, max_cell_types[i]]
        
        cell_type_storage_array.append(max_cell_type)
        cell_value_storage_array.append(max_cell_value)


    cell_type_storage_arrays.append(cell_type_storage_array)    
    cell_value_storage_arrays.append(cell_value_storage_array)     

# Assign to new columns in the dataframe
for i in range(len(cell_type_storage_arrays)):
    merged_df[''.join(['Deconv_cell', str(i + 1)])] = cell_type_storage_arrays[i]
    merged_df[''.join(['Deconv_cell', str(i + 1), '_value'])] = cell_value_storage_arrays[i]


# Since we only consider the top N cell types, we need to correct the weight
# values so that the scatterpies account to the totality of the circle (sum of weights == 1)
deconv_weight_columns = ["Deconv_cell1_value",
                "Deconv_cell2_value",
                "Deconv_cell3_value",
                "Deconv_cell4_value", 
                "Deconv_cell5_value"]

# Create new normalized columns
for i in range(len(cell_type_storage_arrays)):
    
    # Calculate the sum of the top cell type weights
    total = merged_df.loc[:, deconv_weight_columns].sum(axis=1)
    
    # Create column with corrected weight values
    merged_df[''.join(['Deconv_cell', str(i + 1), '_norm_value'])] =  merged_df[''.join(['Deconv_cell', str(i + 1), '_value'])] / total
    

# SLim down the df by selecting columns of interest only
columns_of_interest = ['x', 'y', 
                       "Deconv_cell1", "Deconv_cell1_norm_value",
                       "Deconv_cell2", "Deconv_cell2_norm_value",
                       "Deconv_cell3", "Deconv_cell3_norm_value",
                       "Deconv_cell4", "Deconv_cell4_norm_value", 
                       "Deconv_cell5", "Deconv_cell5_norm_value",]
reduced_df = merged_df.loc[:, columns_of_interest]
        

# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>

# Define color dictionary
colordict = {
    "AC.like": "#CCCCCC",
    "AC.like.Prolif": "#FF6600",
    "Astrocyte": "#00FFCC",
    "B.cell": "#F0E442",
    "CD4.INF": "#0066FF",
    "CD4.rest": "#FF00FF",
    "CD8.cytotoxic": "#00FF00",
    "CD8.EM": "#FF6666",
    "CD8.NK.sig": "#FFCC00",
    "cDC1": "#00FFFF",
    "cDC2": "#FF0066",
    "DC1": "#CCFF00",
    "DC2": "#0000FF",
    "DC3": "#FFCCCC",
    "Endo.arterial": "#CC00FF",
    "Endo.capilar": "#66FF00",
    "Mast": "#FF00CC",
    "MES.like.hypoxia.independent": "#00CCFF",
    "MES.like.hypoxia.MHC": "#003399",
    "Mono.anti.infl": "#FF3366",
    "Mono.hypoxia": "#00FF66",
    "Mono.naive": "#FF9999",
    "Neuron": "#6600FF",
    "NK": "#FFE6E6",
    "NPC.like.neural": "#0072B2",
    "NPC.like.OPC": "#FF0000",
    "NPC.like.Prolif": "#999900",
    "Oligodendrocyte": "#666666",
    "OPC": "#CCFF99",
    "OPC.like": "#000000",
    "OPC.like.Prolif": "#990000",
    "pDC": "#993300",
    "Pericyte": "#996600",
    "Perivascular.fibroblast": "#999999",
    "Plasma.B": "#669900",
    "Prolif.T": "#339900",
    "Reg.T": "#CC79A7",
    "RG": "#009933",
    "Scavenging.endothelial": "#990099",
    "Scavenging.pericyte": "#009900",
    "SMC": "#330099",
    "SMC.COL": "#CC9999",
    "SMC.prolif": "#009999",
    "Stress.sig": "#990066",
    "TAM.BDM.anti.infl": "#990033",
    "TAM.BDM.hypoxia.MES": "#CC3333",
    "TAM.BDM.INF": "#CC6666",
    "TAM.BDM.MHC": "#660099",
    "TAM.MG.aging.sig": "#CCCC99",
    "TAM.MG.pro.infl.I": "#56B4E9",
    "TAM.MG.pro.infl.II": "#333333",
    "TAM.MG.prolif": "#99CC99",
    "Tip.like": "#99CC66",
    "VLMC": "#99CC33"
}



# Smaller df for testing
test_df = reduced_df.iloc[1:1000, ]

# Create dictionary of data to show
data = {
    'x': [ y / 100 for y in test_df.y.tolist()],
    'y': [-x / 100 for x in test_df.x.tolist()],
    'x_full': test_df.y.tolist(),
    'y_full': [ -x for x in test_df.x.tolist()],
    'DeconvCell1': test_df.Deconv_cell1.tolist(), 
    'DeconvCell2': test_df.Deconv_cell2.tolist(), 
    'DeconvCell3': test_df.Deconv_cell3.tolist(), 
    'DeconvCell4': test_df.Deconv_cell4.tolist(), 
    'DeconvCell5': test_df.Deconv_cell5.tolist(),
    'DeconvCell1_w': test_df.Deconv_cell1_norm_value.tolist(),
    'DeconvCell2_w': test_df.Deconv_cell2_norm_value.tolist(),
    'DeconvCell3_w': test_df.Deconv_cell3_norm_value.tolist(),
    'DeconvCell4_w': test_df.Deconv_cell4_norm_value.tolist(),
    'DeconvCell5_w': test_df.Deconv_cell5_norm_value.tolist(),
}

# Convert dictionary to dataframe
df = pd.DataFrame(data)

# Convert dataframe to a ColumnDataSource
source = ColumnDataSource(df)

# Initialize the Bokeh plot
p = figure(width = 1000, height = 1000,
            title = "Deconvolution results", 
            x_axis_label = 'x',
            y_axis_label = 'y', 
            )

for index, row in df.iterrows():
    
    # Get x, y coordinates and category values for the current data point
    x, y = row['x'], row['y']
    categories = row[['DeconvCell1_w', 'DeconvCell2_w', 'DeconvCell3_w', 'DeconvCell4_w', 'DeconvCell5_w']].values
    
    # Extract colors correspoding to cell types
    cell_types = row[['DeconvCell1', 'DeconvCell2', 'DeconvCell3', 'DeconvCell4', 'DeconvCell5']].values
    colors = tuple([colordict[x] for x in cell_types])
    #print(colors)
    
    # Calculate start and end angles for each category
    start_angle = 0
    for i, category_value in enumerate(categories):
        end_angle = start_angle + category_value * 2 * pi
        p.wedge(x=x, y=y, radius=0.05, 
                start_angle=start_angle, end_angle=end_angle,
                line_color="white", fill_color = colors[i], legend_label=f"DeconvCell {i+1}", source=source)
        start_angle = end_angle

# Show no legend
p.legend.visible=False 

# Define tooltip        
TOOLTIPS = [
        ("(x, y)", "(@x_full, @y_full)"), # We show Visium coordinates
        ("TopCell1", "@DeconvCell1, weight = @DeconvCell1_w{0.00}"),
        ("TopCell2", "@DeconvCell2, weight = @DeconvCell2_w{0.00}"),
        ("TopCell3", "@DeconvCell3, weight = @DeconvCell3_w{0.00}"),
        ("TopCell4", "@DeconvCell4, weight = @DeconvCell4_w{0.00}"),
        ("TopCell5", "@DeconvCell5, weight = @DeconvCell5_w{0.00}"),        

]

# Assign tooltips to HoverTool
hover = HoverTool(tooltips = TOOLTIPS)        

# Add the hover tool to the plot
p.add_tools(hover)

# Output to HTML file
output_file("visium_deconvolution_plot_test.html")

# Save the plot
save(p)




