# -*- coding: utf-8 -*-
"""

Script to generate an interactive plot to show the top n cell types obtained
from an spatial deconvolution analysis
See documenation for use parameters .
"""
# <! ------------------------------------------------------------------------!>
# <!                       IMPORTS                                           !>
# <! ------------------------------------------------------------------------!>
import pandas as pd
from math import pi
from bokeh.plotting import figure, output_file, save,output_notebook, show, curdoc

from bokeh.transform import cumsum
from bokeh.models import ColumnDataSource, HoverTool,Range1d
from bokeh.palettes import Category20
from bokeh.models import CustomJS, TapTool
from PIL import Image
import base64
import io
# Define color dictionary
clusters_colordict = {
    0: "#CCCCCC",
    1: "#FF6600",
    2: "#00FFCC",
    3: "#F0E442",
    4: "#0066FF",
    5: "#FF00FF",
    6: "#00FF00",
    7: "#FF6666",
    8: "#FFCC00",
    9: "#00FFFF",
    10: "#FF0066",
    11: "#CCFF00",
    12: "#0000FF",
    13: "#FFCCCC",
    14: "#CC00FF",
}

# Define color dictionary
colordict = {
    "AClike": "#CCCCCC",
    "AClikeProlif": "#FF6600",
    "Astrocyte": "#00FFCC",
    "Bcell": "#F0E442",
    "CD4INF": "#0066FF",
    "CD4rest": "#FF00FF",
    "CD8cytotoxic": "#00FF00",
    "CD8EM": "#FF6666",
    "CD8NKsig": "#FFCC00",
    "cDC1": "#00FFFF",
    "cDC2": "#FF0066",
    "DC1": "#CCFF00",
    "DC2": "#0000FF",
    "DC3": "#FFCCCC",
    "Endoarterial": "#CC00FF",
    "Endocapilar": "#66FF00",
    "Mast": "#FF00CC",
    "MESlikehypoxiaindependent": "#00CCFF",
    "MESlikehypoxiaMHC": "#003399",
    "Monoantiinfl": "#FF3366",
    "Monohypoxia": "#00FF66",
    "Mononaive": "#FF9999",
    "Neuron": "#6600FF",
    "NK": "#FFE6E6",
    "NPClikeneural": "#0072B2",
    "NPClikeOPC": "#FF0000",
    "NPClikeProlif": "#999900",
    "Oligodendrocyte": "#666666",
    "OPC": "#CCFF99",
    "OPClike": "#000000",
    "OPClikeProlif": "#990000",
    "pDC": "#993300",
    "Pericyte": "#996600",
    "Perivascularfibroblast": "#999999",
    "PlasmaB": "#669900",
    "ProlifT": "#339900",
    "RegT": "#CC79A7",
    "RG": "#009933",
    "Scavengingendothelial": "#990099",
    "Scavengingpericyte": "#009900",
    "SMC": "#330099",
    "SMCCOL": "#CC9999",
    "SMCprolif": "#009999",
    "Stresssig": "#990066",
    "TAMBDMantiinfl": "#990033",
    "TAMBDMhypoxiaMES": "#CC3333",
    "TAMBDMINF": "#CC6666",
    "TAMBDMMHC": "#660099",
    "TAMMGagingsig": "#CCCC99",
    "TAMMGproinflI": "#56B4E9",
    "TAMMGproinflII": "#333333",
    "TAMMGprolif": "#99CC99",
    "Tiplike": "#99CC66",
    "VLMC": "#99CC33",
    "MESlikehypoxiaindependent" : "#990033",
}

# <! ------------------------------------------------------------------------!>
# <!                           DATA PREPARATION                              !>
# <! ------------------------------------------------------------------------!>

def process_data(norm_weights_filepath, st_coords_filepath, data_clustered, deconv_method, n_largest_cell_types, scale_factor):
    # Read spatial deconvolution result CSV file
    norm_weights_df = pd.read_csv(norm_weights_filepath, sep = '\t')
    norm_weights_df.index.name = None

    # Read spatial coordinates CSV file
    st_coords_df = pd.read_csv(st_coords_filepath, header=None).set_index(0)
    st_coords_df.index.name = None
    st_coords_df.columns = [ "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    st_coords_df["pxl_row_in_fullres"] = st_coords_df["pxl_row_in_fullres"]*scale_factor
    st_coords_df["pxl_col_in_fullres"] = st_coords_df["pxl_col_in_fullres"]*scale_factor
    # It will be difficult to show the information of all 54 cell types when hovering
    # Thus, for each barcoded spot, retrieve the maximum 5 weights and create new columns
    # accordingly. Those 5 max columns will be the info shown in the hovertool
    max_weights = norm_weights_df.apply(lambda x: x.nlargest(n_largest_cell_types).index.values, axis=1)

    merged_df = pd.concat([st_coords_df, norm_weights_df], axis = 1, join = 'inner')

    data_with_clusters = pd.read_csv(data_clustered)
    clusters_col =  pd.DataFrame(data_with_clusters["BayesSpace"]).set_index(data_with_clusters["Unnamed: 0"])
    merged_df["Cluster"] = clusters_col


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

    # print(len(cell_type_storage_arrays[0]))
    # Assign to new columns in the dataframe
    for i in range(n_largest_cell_types):
        merged_df[''.join([deconv_method,  '_Deconv_cell', str(i + 1)])] = cell_type_storage_arrays[i]
        merged_df[''.join([deconv_method, '_Deconv_cell', str(i + 1), '_value'])] = cell_value_storage_arrays[i]

    # Since we only consider the top N cell types, we need to correct the weight
    # values so that the scatterpies account to the totality of the circle (sum of weights == 1)

    deconv_weight_columns = [f"{deconv_method}_Deconv_cell{i + 1}_value" for i in range(n_largest_cell_types)]

    # Create new normalized columns
    for i in range(n_largest_cell_types):

        # Calculate the sum of the top cell type weights
        total = merged_df.loc[:, deconv_weight_columns].sum(axis=1)

        # Create column with corrected weight values
        merged_df[''.join([deconv_method, '_Deconv_cell', str(i + 1), '_norm_value'])] =  merged_df[''.join([deconv_method, '_Deconv_cell', str(i + 1), '_value'])] / total

    # SLim down the df by selecting columns of interest only
    columns_of_interest = ['pxl_row_in_fullres', 'pxl_col_in_fullres','Cluster' , "in_tissue"] + [f"{deconv_method}_Deconv_cell{i + 1}_norm_value" for i in range(n_largest_cell_types)] \
        + [f"{deconv_method}_Deconv_cell{i + 1}" for i in range(n_largest_cell_types)]
    reduced_df = merged_df.loc[:, columns_of_interest]
    return reduced_df


import base64
import io
import math

def image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode()
    return f"data:image/png;base64,{encoded_string}"
def get_image_display_infos(image_path):
    im = Image.open(image_path).convert("RGB")
    # Merge coordinate df and cell weight df
    image_display_infos = {
        "image_path" : image_path,
        "x0" : 0,
        "y0" : 0,
        "im_w" : im.size[0],
        "im_h" : im.size[1],

    }
    return image_display_infos


def calculate_stddev(vector):
    mean = sum(vector) / len(vector)
    # Calculate the variance
    variance = sum((x - mean) ** 2 for x in vector) / len(vector)
    stddev = math.sqrt(variance)
    return stddev

def calculate_rmsd(vector1, vector2):
    if len(vector1) != len(vector2):
        raise ValueError("The vectors must be of the same length")
    
    squared_diffs = [(v1 - v2) ** 2 for v1, v2 in zip(vector1, vector2)]
    mean_squared_diff = sum(squared_diffs) / len(vector1)
    return math.sqrt(mean_squared_diff)


def post_process_data(norm_weights_filepaths, st_coords_filepath, data_clustered, deconv_methods, n_largest_cell_types, scale_factor):
    norm_weights_dfs = [process_data(props, st_coords_filepath,data_clustered, deconv_methods[index], n_largest_cell_types, scale_factor = scale_factor)\
                      for index, props  in enumerate(norm_weights_filepaths)]
    processed_data = norm_weights_dfs[0]
    for i in range(1, len(norm_weights_dfs)):
        processed_data = pd.merge(processed_data, norm_weights_dfs[i])
    nb_spots_samples = processed_data.shape[0]

    for method in deconv_methods:
        processed_data[f'{method}_values'] = processed_data.apply(
            lambda row: [row[f"{method}_Deconv_cell{i}_norm_value"] for i in range(1, n_largest_cell_types + 1)],
            axis=1
        )
    #process data again because i dont want to change all process_data
    # i need methods comparing data, so i had to add this step

    weights_gathered = pd.read_csv(norm_weights_filepaths[0], sep = '\t')
    weights_gathered.index.name = None

    cell_types = weights_gathered.columns.copy()
    weights_gathered  =  weights_gathered.rename(columns=lambda x: f'{deconv_methods[0]}_{x}')
    for i in range(1,len(norm_weights_filepaths)):
      tmp_df = pd.read_csv(norm_weights_filepaths[i], sep = '\t')
      tmp_df.index.name = None
      tmp_df  =  tmp_df.rename(columns=lambda x: f'{deconv_methods[i]}_{x}')
      weights_gathered = pd.merge(weights_gathered, tmp_df, left_index=True, right_index=True)
    def comparing_value_stdev(row):
      st_devs = []
      for t in cell_types:
        st_dev = calculate_stddev( [ row[f"{method}_{t}"] for method in deconv_methods] )
        st_devs.append(st_dev)
      return calculate_stddev(st_devs)


    weights_gathered.reset_index(drop=True, inplace = True)
    def comparing_value_rmsd(row):
      v0 = [  row[f"{deconv_methods[0]}_{t}"] for t in cell_types]
      v1 = [  row[f"{deconv_methods[1]}_{t}"] for t in cell_types]
      return calculate_rmsd(v0 , v1)

    if len(deconv_methods) == 2:
      weights_gathered["error_value"] = weights_gathered.apply(lambda row : comparing_value_rmsd(row), axis = 1)
    elif len(deconv_methods) > 2:
      weights_gathered["error_value"] = weights_gathered.apply(lambda row : comparing_value_stdev(row), axis = 1)
    processed_data["error_value"] = weights_gathered["error_value"]
    return processed_data

# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>
from bokeh.events import ButtonClick
from bokeh.models import BoxAnnotation, Label, Plot, Rect, Text, Button, CustomJS, Div,Slider, PanTool
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from bokeh.layouts import column, row, gridplot,Spacer
from PIL import Image
import numpy as np

import math

    
   
def vis_with_separate_clusters_view(reduced_df, image_path, deconv_methods, nb_spots_samples, output , show_legend = False, show_figure = False ):
    image_display_infos = get_image_display_infos(image_path)
    image_display_infos = {x: int(np.ceil(image_display_infos[x]/2)) if x != "image_path" else image_display_infos[x] for x in image_display_infos}
    # Smaller sample
    test_df = reduced_df[reduced_df["in_tissue"] == 1].head(nb_spots_samples).copy()
    # Create a single tooltip column for each circle
    test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'.join( \
                                            [f"<span style='color: red;'> Spot</span> : (x = { row['pxl_col_in_fullres']/2:.2f}, y = {-row['pxl_row_in_fullres']/2:.2f})"] ),\
                                            axis=1)
    test_df['error_tooltip_data'] = test_df.apply(lambda row: '<br>'.join( \
                                            [f"<span style='color: red;'> Spot</span> : (x = { row['pxl_col_in_fullres']/2:.2f}, y = {-row['pxl_row_in_fullres']/2:.2f})"]\
                                                + [f"<span style='color: blue;'> Cluster</span> : {row['Cluster']}"]),\
                                            axis=1)
    # Update the data dictionary
    data = {
        'x': [y/2 for y in test_df.pxl_col_in_fullres.tolist()],
        'y': [-x/2 for x in test_df.pxl_row_in_fullres.tolist()],
        'tooltip_data': test_df['tooltip_data'].tolist(),
        'Cluster' : test_df['Cluster'].tolist() ,
        'error_value' : test_df["error_value"].tolist(),
        'error_tooltip_data': test_df['error_tooltip_data'].tolist()
    }
    # Convert dictionary to dataframe
    df = pd.DataFrame(data)
    # Initialize the Bokeh plot
    p = figure(width = image_display_infos.get("im_w"), height = image_display_infos.get("im_h"),
                title = "Clustering results",
                x_axis_label = 'x',
                y_axis_label = 'y',
                output_backend="webgl"
                )
    # Add the image with a ColumnDataSource
    image_source = ColumnDataSource(data=dict(
        url=[ image_to_base64(image_display_infos.get("image_path"))],
        x=[ image_display_infos.get("x0") ],
        y=[ image_display_infos.get("y0") ],
        w=[image_display_infos.get("im_w") ],
        h=[ image_display_infos.get("im_h") ],
        alpha=[1.0]  # Initial alpha value
    ))
    image = p.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)
    # Create a slider for image transparency
    slider = Slider(start=0, end=1, value=1, step=.1, title="Image Transparency")
    # Create a callback to update the image alpha
    callback = CustomJS(args=dict(image_source=image_source), code="""
        var alpha = cb_obj.value;
        image_source.data['alpha'] = [alpha];
        image_source.change.emit();
    """)

    slider.js_on_change('value', callback)
    # Create a dictionary to store scatter renderers
    scatter_renderers = {}
    # Group the dataframe by cluster
    grouped = df.groupby('Cluster')
    # Plot each cluster separately
    for cluster, group in grouped:
        color = clusters_colordict.get(cluster, '#000000')
        source = ColumnDataSource(group)

        scatter = p.scatter(x='x', y='y', size=5,
                            marker="circle",  # Specify the marker shape
                            fill_color=color
                            , line_width=0,
                            source=source,
                            legend_label=f"Cluster {int(cluster)}")
        scatter_renderers[cluster] = scatter
    for method in deconv_methods:
        # Create a single tooltip column for each circle
        test_df[f"{method}_tooltip_data"] = test_df.apply(lambda row: '<br>'.join([
            f"<div style='display:flex;align-items:center;'>"
            f"<div style='width:10px;height:10px;background-color:{colordict.get(row[f'{method}_Deconv_cell{i+1}'], '#000000')};margin-right:5px;'></div>"
            f"<span style='color: blue;'>{row[f'{method}_Deconv_cell{i+1}']}</span>: {row[f'{method}_Deconv_cell{i+1}_norm_value']*100:.2f}%"
            f"</div>"
            for i in range(n_largest_cell_types)
        ] +  [f"<span style='color: red;'> Spot</span> : (x = {row['pxl_col_in_fullres']/2:.2f}, y = {-row['pxl_row_in_fullres']/2:.2f})"]), axis=1)
        data[f"{method}_tooltip_data"] = test_df[f"{method}_tooltip_data"].tolist()
        for i in range(1, n_largest_cell_types + 1):
            data[f'{method}_DeconvCell{i}'] = test_df[f'{method}_Deconv_cell{i}'].tolist()
            data[f'{method}_DeconvCell{i}_w'] = test_df[f'{method}_Deconv_cell{i}_norm_value'].tolist()
    # Convert dictionary to dataframe
    df = pd.DataFrame(data)
    deconv_plots = []
    # print(df.head(5))
    for method in deconv_methods:
        plot = figure(width=image_display_infos.get("im_w"),
                    height=image_display_infos.get("im_h"),
                title="Deconvolution results",
                x_axis_label='x',
                y_axis_label='y',
                output_backend="webgl",
                )
        plot.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)
        # Create a Div for displaying the message
        for index, row in df.iterrows():
            x, y = row['x'], row['y']
            categories = row[[f'{method}_DeconvCell{i+1}_w' for i in range(n_largest_cell_types)]].values
            cell_types = row[[f'{method}_DeconvCell{i+1}' for i in range(n_largest_cell_types)]].values
            colors = tuple([colordict[x] for x in cell_types])
            # Create a single ColumnDataSource for all wedges in this circle
            circle_source = ColumnDataSource({
                'x': [x],
                'y': [y],
                f"{method}_tooltip_data": [row[f"{method}_tooltip_data"]]
            })
            start_angle = 0
            for i, category_value in enumerate(categories):
                end_angle = start_angle + category_value * 2 * pi
                wedge = plot.wedge(x='x', y='y', radius=4.7,
                        start_angle=start_angle, end_angle=end_angle,
                        line_width=0, fill_color=colors[i],
                        legend_label=f"Cluster {row['Cluster']}", source=circle_source, visible=False)
                start_angle = end_angle
        deconv_plots.append(plot)
    from bokeh.models import LinearColorMapper, ColorBar
    from bokeh.transform import transform
    from bokeh.palettes import Viridis256

    rmsd_plot = figure(width=image_display_infos.get("im_w"),
                        height=image_display_infos.get("im_h"),
                        title="Deconvolution results comparing",
                        x_axis_label='x',
                        y_axis_label='y',
                        output_backend="webgl")

    rmsd_plot.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)
    min_val, max_val = min(data["error_value"]), max(data["error_value"])
    color_map = LinearColorMapper(palette=Viridis256, low=min_val, high=max_val)
    for index, row in df.iterrows():
        x, y = row['x'], row['y']
        error_value = row['error_value']
        error_tooltip_data = row["error_tooltip_data"]
        circle_source = ColumnDataSource({
            'x': [x],
            'y': [y],
            'error_value': [error_value] , # Add rmsd_value for color mapping
            'error_tooltip_data' : [error_tooltip_data]
        })

        rmsd_plot.scatter(x='x', y='y', size=5,
                            marker="circle",
                            fill_color=transform('error_value', color_map),  # Correct color mapping
                            line_width=0,
                            source=circle_source,
                            visible=True)

    # Add a color bar for the color mapping
    color_bar = ColorBar(color_mapper=color_map,
                        label_standoff=14,
                        location=(0, 0),
                        title='Color Range')

    # Add color bar to the plot
    rmsd_plot.add_layout(color_bar, 'right')
    ########

    text1 = """
        <div style="
            background-color: #f0f0f0;
            border: 2px solid #3c3c3c;
            border-radius: 10px;
            padding: 15px;
            margin: 10px 0;
            font-family: 'Helvetica', 'Arial', sans-serif;
            font-size: 14px;
            line-height: 1.5;
            color: #333333;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        ">
            <h3 style="
                margin-top: 0;
                color: #2c3e50;
                font-size: 18px;
                border-bottom: 1px solid #bdc3c7;
                padding-bottom: 8px;
            ">Visualization Infos</h3>
            <p style="margin-bottom: 0;">
                This view shows the clusters. Each point represents a spot,
                    and the colors indicate the different clusters."
            </p>
        </div>
        """
    text2 = """
        <div style="
            background-color: #f0f0f0;
            border: 2px solid #3c3c3c;
            border-radius: 10px;
            padding: 15px;
            margin: 10px 0;
            font-family: 'Helvetica', 'Arial', sans-serif;
            font-size: 14px;
            line-height: 1.5;
            color: #333333;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        ">
            <h3 style="
                margin-top: 0;
                color: #2c3e50;
                font-size: 18px;
                border-bottom: 1px solid #bdc3c7;
                padding-bottom: 8px;
            ">Visualization Infos</h3>
            <p style="margin-bottom: 0;">
                This visualization shows deconvolution results by cluster. 
                In each spot, celltype proportions are represented with a pie chart"
            </p>
        </div>
        """
    text3 = """
        <div style="
            background-color: #f0f0f0;
            border: 2px solid #3c3c3c;
            border-radius: 10px;
            padding: 15px;
            margin: 10px 0;
            font-family: 'Helvetica', 'Arial', sans-serif;
            font-size: 14px;
            line-height: 1.5;
            color: #333333;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        ">
            <h3 style="
                margin-top: 0;
                color: #2c3e50;
                font-size: 18px;
                border-bottom: 1px solid #bdc3c7;
                padding-bottom: 8px;
            ">Visualization Infos</h3>
            <p style="margin-bottom: 0;">
                This visualization compares deconvolution results between different methods. 
                When having 2 methods the comparing is based on RMSD metric. However, when 
                using more than 2 methods, for each spot, and for each celltype, standard 
                deviation of deconvolution values among methods is calculated, then the final
                error value for that spot is the standard deviation of calculated standard
                    deviations among celltypes"
            </p>
        </div>
        """
    title_text = """
        <div style="
            background-color: #ffffff;
            border: 3px solid #2c3e50;
            border-radius: 15px;
            padding: 30px;
            margin: 20px 0;
            font-family: 'Helvetica', 'Arial', sans-serif;
            color: #333333;
            text-align: center;
            box-shadow: 0 8px 16px rgba(0, 0, 0, 0.1);
        ">
            <h1 style="
                margin-top: 0;
                font-size: 25px;
                color: #2c3e50;
                letter-spacing: 1px;
                font-weight: bold;
            ">Spatial Transcriptomic Exploration Interface</h1>
            <hr style="
                width: 80px;
                border: 2px solid #2c3e50;
                margin: 10px auto;
            ">
        </div>
        """
    # Créez le widget Div
    info_box = Div(
        text= text1,
        width=700,
        height=200
    )

    # Modifiez les callbacks des boutons pour mettre à jour le texte du Div
    show_all_button = Button(label="Show Clusters", width=100, button_type = 'primary')
    show_all_button.js_on_click(CustomJS(args=dict(p=p,  rmsd_plot = rmsd_plot, deconv_plots=deconv_plots, info_box=info_box, text1 = text1), code="""
        p.visible = true;
        deconv_plots.forEach((p) => {p.visible = false});
        rmsd_plot.visible = false;
        info_box.text = text1;
    """))
    spacer = Spacer(width=50)  # Adjust the width as needed
    button_methods = []
    for index, method in enumerate(deconv_methods):
        button = Button(label=f"Show deconvolution by Cluster with {method}", width=150, button_type = 'primary')
        button.js_on_click(CustomJS(args=dict(p=p,  rmsd_plot = rmsd_plot,  deconv_plots=deconv_plots, index = index, info_box=info_box, text2 = text2), code="""
            p.visible = false;
            rmsd_plot.visible = false;
            deconv_plots.forEach((p) => {p.visible = false});
            deconv_plots[index].visible = true;
            info_box.text = text2;
        """))
        button_methods.append(button)
        spacer1 = Spacer(width=150)  # Adjust the width as needed
        button_methods.append(spacer1)

    rmsd_button = Button(label="Compare Deconvolution Methods", width=100, button_type = 'primary')
    rmsd_button.js_on_click(CustomJS(args=dict(p=p, rmsd_plot = rmsd_plot, deconv_plots=deconv_plots, info_box=info_box, text3 = text3), code="""
        p.visible = false;
        deconv_plots.forEach((p) => {p.visible = false});
        rmsd_plot.visible = true; 
        info_box.text = text3;

    """))

    spacer1 = Spacer(width=100)
    spacer2 = Spacer(width=100)
    spacer3 = Spacer(width=100)
    spacer4 = Spacer(width=100)

    # Assuming you have your data in a pandas DataFrame called 'df'
    csv_source = ColumnDataSource({'data': [df.drop(columns=[f"{method}_tooltip_data" for method in deconv_methods]).to_csv(index=False)]})
    download_button = Button(label="Download raw data", width=100, button_type = 'primary')
    download_button.js_on_click(CustomJS(args=dict(source=csv_source), code="""
        const data = source.data['data'][0];
        const blob = new Blob([data], { type: 'text/csv;charset=utf-8;' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = "raw_data.csv";
        link.click();
        URL.revokeObjectURL(url);
    """))
    tap_tool = TapTool()
    p.add_tools(tap_tool)
    # Associez le callback au TapTool
    hover = HoverTool(tooltips="""
        <div style="width:220px">
            @tooltip_data
        </div>
    """)
    hover_error = HoverTool(tooltips="""
        <div style="width:220px">
            @error_tooltip_data
        </div>
    """)
    rmsd_plot.add_tools(hover_error)
    hover_list = []
    for method in deconv_methods:
        tooltip_string = f"""
            <div style="width:220px">
                @{method}_tooltip_data
            </div>
        """
        # Create the HoverTool with the constructed tooltip string
        h = HoverTool(tooltips=tooltip_string)
        hover_list.append(h)
    p.add_tools(hover)
    for index, plot in enumerate(deconv_plots):
        plot.add_tools(hover_list[index])
        leg_plot = plot.legend[0]
        leg_plot.glyph_width = 0
        leg_plot.label_text_font_size = "15pt"
        plot.add_layout(leg_plot,'right')
        plot.legend.location = "top_right"
        plot.legend.click_policy = "hide"
        plot.visible = False
        plot.legend.visible = True
        sorted_items_plot = sorted(leg_plot.items, key=lambda item: item.label['value'])
        leg_plot.items = sorted_items_plot
    leg = p.legend[0]
    leg.glyph_height = 20
    leg.glyph_width = 20
    leg.label_text_font_size = "15pt"
    sorted_items = sorted(leg.items, key=lambda item: item.label['value'])
    leg.items = sorted_items
    p.add_layout(leg,'right')

    title_box = Div(
        text= title_text,
        width=700,
        height=100
    )

    from bokeh.layouts import column, row

    # Set size for the main plots
    p.width = 900  # Set the width
    p.height = 700  # Set the height
    rmsd_plot.width = 900  # Set the width
    rmsd_plot.height = 700  # Set the height

    # Set size for each plot in 'deconv_plots'
    for plot in deconv_plots:
        plot.width = 900 
        plot.height = 700  

    # Create the row for the buttons and slider
    buttons_col = column(
        show_all_button, *button_methods,
        rmsd_button, download_button , slider,
        Spacer(height=20)  
    )

    controls_column = column(
        title_box,
        Spacer(height = 60),
        info_box,
        buttons_col,
    )
    plots_column = column(
        p,            
        *deconv_plots, 
        rmsd_plot,
    )
    layout = row(
        controls_column,  
        plots_column,    
        sizing_mode='stretch_both'  # Stretch both components to fit the available space
    )

    p.legend.location = "top_right"
    p.legend.click_policy = "hide"
    p.legend.visible = True
    rmsd_plot.visible = False

    if show_figure:
        show(layout)
    output_file(output, mode='inline')
    save(layout)


# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>

if __name__ == "__main__":
    import sys
    argv = sys.argv
    sp_input = argv[1]
    norm_weights_filepaths = argv[2].split(',')
    st_coords_filepath = argv[3]
    data_clustered = argv[4]
    image_path = argv[5]
    n_largest_cell_types = int(argv[6])
    scale_factor = float(argv[7])
    output_html = argv[8]
    deconv_methods = argv[9].split(',')
    print("Processing data ...\n")
    processed_data = post_process_data(norm_weights_filepaths=norm_weights_filepaths, st_coords_filepath=st_coords_filepath, data_clustered=data_clustered, \
                                 deconv_methods=deconv_methods, n_largest_cell_types=n_largest_cell_types, scale_factor=scale_factor)

    print(f"Deconvolution methods {deconv_methods}")
    nb_spots_samples = processed_data.shape[0]
    print(f"Generating vis with {nb_spots_samples} spots and top {n_largest_cell_types} cells...\n")
    vis_with_separate_clusters_view(reduced_df=processed_data, image_path = image_path, deconv_methods = deconv_methods,nb_spots_samples = nb_spots_samples, output= output_html )
