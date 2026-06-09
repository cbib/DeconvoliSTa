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
from bokeh.models import BoxAnnotation, Label, Plot, Rect, Text, Button, CustomJS, Div,Slider, PanTool, CheckboxButtonGroup
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from bokeh.layouts import column, row, gridplot,Spacer
from PIL import Image
import numpy as np

import math

    
   
def vis_with_separate_clusters_view(reduced_df, image_path, deconv_methods, nb_spots_samples, n_largest_cell_types, output, cluster_composition=None, show_legend=False, show_figure=False):
    from bokeh.models import LinearColorMapper, ColorBar
    from bokeh.palettes import Viridis256

    image_display_infos = get_image_display_infos(image_path)
    image_display_infos = {x: int(np.ceil(image_display_infos[x]/2)) if x != "image_path" else image_display_infos[x] for x in image_display_infos}

    test_df = reduced_df[reduced_df["in_tissue"] == 1].head(nb_spots_samples).copy().reset_index(drop=True)
    n_spots = len(test_df)

    all_x = [y/2 for y in test_df.pxl_col_in_fullres.tolist()]
    all_y = [-x/2 for x in test_df.pxl_row_in_fullres.tolist()]

    test_df['tooltip_data'] = test_df.apply(
        lambda row: f"<span style='color: red;'> Spot</span> : (x = {row['pxl_col_in_fullres']/2:.2f}, y = {-row['pxl_row_in_fullres']/2:.2f})",
        axis=1)
    test_df['error_tooltip_data'] = test_df.apply(
        lambda row: (f"<span style='color: red;'> Spot</span> : (x = {row['pxl_col_in_fullres']/2:.2f}, y = {-row['pxl_row_in_fullres']/2:.2f})<br>"
                     f"<span style='color: blue;'> Cluster</span> : {row['Cluster']}"),
        axis=1)

    # Shared image source (one object, reused across all plots)
    image_source = ColumnDataSource(data=dict(
        url=[image_to_base64(image_display_infos.get("image_path"))],
        x=[image_display_infos.get("x0")],
        y=[image_display_infos.get("y0")],
        w=[image_display_infos.get("im_w")],
        h=[image_display_infos.get("im_h")],
        alpha=[1.0]
    ))

    slider = Slider(start=0, end=1, value=1, step=.1, title="Image Transparency")
    slider.js_on_change('value', CustomJS(args=dict(image_source=image_source), code="""
        var alpha = cb_obj.value;
        image_source.data['alpha'] = [alpha];
        image_source.change.emit();
    """))

    # --- Cluster plot (one scatter renderer per cluster, already O(n_clusters)) ---
    p = figure(width=900, height=700, title="Clustering results",
               x_axis_label='x', y_axis_label='y')
    p.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)

    cluster_source_df = pd.DataFrame({
        'x': all_x, 'y': all_y,
        'Cluster': test_df['Cluster'].tolist(),
        'tooltip_data': test_df['tooltip_data'].tolist()
    })
    cluster_renderers = []
    cluster_ids = []
    for cluster, group in cluster_source_df.groupby('Cluster'):
        color = clusters_colordict.get(cluster, '#000000')
        r = p.scatter(x='x', y='y', size=5, marker="circle", fill_color=color,
                  line_width=0, source=ColumnDataSource(group),
                  legend_label=f"Cluster {int(cluster)}")
        cluster_renderers.append(r)
        cluster_ids.append(int(cluster))
    p.add_tools(HoverTool(tooltips="<div style='width:220px'>@tooltip_data{safe}</div>"))

    # --- Deconv plots: vectorized — 1 source per method, n_largest_cell_types renderers ---
    # Replaces the previous n_spots × n_cell_types individual sources/renderers loop.
    deconv_plots = []
    deconv_sources = []
    for method in deconv_methods:
        test_df[f"{method}_tooltip_data"] = test_df.apply(lambda row: '<br>'.join([
            f"<div style='display:flex;align-items:center;'>"
            f"<div style='width:10px;height:10px;background-color:{colordict.get(row[f'{method}_Deconv_cell{i+1}'], '#000000')};margin-right:5px;'></div>"
            f"<span style='color: blue;'>{row[f'{method}_Deconv_cell{i+1}']}</span>: {row[f'{method}_Deconv_cell{i+1}_norm_value']*100:.2f}%"
            f"</div>"
            for i in range(n_largest_cell_types)
        ] + [f"<span style='color: red;'> Spot</span> : (x = {row['pxl_col_in_fullres']/2:.2f}, y = {-row['pxl_row_in_fullres']/2:.2f})"]), axis=1)

        # Cumulative angles per spot: shape (n_spots, n_largest_cell_types+1)
        weights = np.array([
            [test_df.iloc[i][f'{method}_Deconv_cell{j+1}_norm_value'] for j in range(n_largest_cell_types)]
            for i in range(n_spots)
        ])
        cumulative = np.hstack([np.zeros((n_spots, 1)), np.cumsum(weights, axis=1)]) * 2 * pi

        # One ColumnDataSource for all spots × all cell type positions
        source_data = {
            'x': all_x,
            'y': all_y,
            'cluster': [int(c) for c in test_df['Cluster'].tolist()],
            'alpha': [1.0] * n_spots,
            'tooltip_data': test_df[f"{method}_tooltip_data"].tolist()
        }
        for j in range(n_largest_cell_types):
            source_data[f'start_{j}'] = cumulative[:, j].tolist()
            source_data[f'end_{j}'] = cumulative[:, j + 1].tolist()
            source_data[f'color_{j}'] = [
                colordict.get(test_df.iloc[i][f'{method}_Deconv_cell{j+1}'], '#000000')
                for i in range(n_spots)
            ]
        shared_source = ColumnDataSource(source_data)

        plot = figure(width=900, height=700, title=f"Deconvolution results - {method}",
                      x_axis_label='x', y_axis_label='y',
                      x_range=p.x_range, y_range=p.y_range)
        plot.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)

        for j in range(n_largest_cell_types):
            plot.wedge(x='x', y='y', radius=4.7,
                       start_angle=f'start_{j}', end_angle=f'end_{j}',
                       fill_color=f'color_{j}', fill_alpha='alpha', line_width=0, source=shared_source)

        plot.add_tools(HoverTool(tooltips="<div style='width:220px'>@tooltip_data{safe}</div>"))
        plot.visible = False
        deconv_plots.append(plot)
        deconv_sources.append(shared_source)

    # --- RMSD plot: single vectorized scatter (replaces n_spots individual sources) ---
    error_values = test_df["error_value"].tolist()
    min_val, max_val = min(error_values), max(error_values)
    color_map = LinearColorMapper(palette=Viridis256, low=min_val, high=max_val)

    rmsd_source = ColumnDataSource({
        'x': all_x, 'y': all_y,
        'cluster': [int(c) for c in test_df['Cluster'].tolist()],
        'alpha': [1.0] * n_spots,
        'error_value': error_values,
        'error_tooltip_data': test_df['error_tooltip_data'].tolist()
    })
    rmsd_plot = figure(width=900, height=700, title="Deconvolution results comparing",
                       x_axis_label='x', y_axis_label='y',
                       x_range=p.x_range, y_range=p.y_range)
    rmsd_plot.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)
    rmsd_plot.scatter(x='x', y='y', size=5, marker="circle",
                      fill_color={'field': 'error_value', 'transform': color_map},
                      fill_alpha='alpha', line_width=0, source=rmsd_source)
    rmsd_plot.add_tools(HoverTool(tooltips="<div style='width:220px'>@error_tooltip_data{safe}</div>"))

    color_bar = ColorBar(color_mapper=color_map, label_standoff=14, location=(0, 0), title='Color Range')
    rmsd_plot.add_layout(color_bar, 'right')
    rmsd_plot.visible = False

    # --- Cluster filter SHARED across all views (clusters + each method + comparison) ---
    # A cluster selection applies everywhere (like the zoom). Select/Deselect all buttons.
    cluster_filter_code = """
        const active = new Set(checkbox.active.map(i => cluster_ids[i]));
        for (let k = 0; k < cluster_renderers.length; k++) {
            cluster_renderers[k].visible = active.has(cluster_ids[k]);
        }
        const all_sources = deconv_sources.concat([rmsd_source]);
        for (const s of all_sources) {
            const cl = s.data['cluster'];
            const al = s.data['alpha'];
            for (let i = 0; i < cl.length; i++) { al[i] = active.has(cl[i]) ? 1 : 0; }
            s.change.emit();
        }
    """
    cluster_checkbox = CheckboxButtonGroup(
        labels=[f"Cluster {c}" for c in cluster_ids],
        active=list(range(len(cluster_ids)))
    )
    cluster_checkbox.js_on_change('active', CustomJS(
        args=dict(checkbox=cluster_checkbox, deconv_sources=deconv_sources,
                  rmsd_source=rmsd_source, cluster_renderers=cluster_renderers,
                  cluster_ids=cluster_ids),
        code=cluster_filter_code))
    select_all_btn = Button(label="Select all clusters", width=120, button_type='success')
    select_all_btn.js_on_click(CustomJS(
        args=dict(checkbox=cluster_checkbox, cluster_ids=cluster_ids),
        code="checkbox.active = Array.from(Array(cluster_ids.length).keys());"))
    deselect_all_btn = Button(label="Deselect all clusters", width=120, button_type='warning')
    deselect_all_btn.js_on_click(CustomJS(
        args=dict(checkbox=cluster_checkbox),
        code="checkbox.active = [];"))
    cluster_controls = column(
        Div(text="<b>Clusters</b> (filter applied to all views):", width=320),
        row(select_all_btn, deselect_all_btn),
        cluster_checkbox
    )

    # --- Average cell-type composition panel for the selected clusters ---
    comp_div = Div(text="<i>Select one or more cluster(s) to see their cell-type composition.</i>",
                   width=560, height=460)
    if cluster_composition is not None:
        comp_cb = CustomJS(args=dict(checkbox=cluster_checkbox, cluster_ids=cluster_ids,
                                     comp=cluster_composition, div=comp_div,
                                     methods=deconv_methods, colordict=colordict),
            code="""
            const active = checkbox.active.map(i => cluster_ids[i]);
            const cts = comp['celltypes'];
            if (active.length === 0) { div.text = "<i>No cluster selected.</i>"; return; }

            function methodBlock(m, clusterList) {
                const means = comp['means'][m];
                const counts = comp['counts'][m];
                let agg = new Array(cts.length).fill(0);
                let tw = 0;
                for (const cl of clusterList) {
                    const key = String(cl);
                    if (!(key in means)) continue;
                    const w = counts[key];
                    const v = means[key];
                    for (let j = 0; j < cts.length; j++) { agg[j] += w * v[j]; }
                    tw += w;
                }
                if (tw === 0) return "";
                for (let j = 0; j < cts.length; j++) { agg[j] /= tw; }
                let idx = Array.from(cts.keys()).sort((a, b) => agg[b] - agg[a]).slice(0, 4);
                let h = "<div style='margin-top:2px'><u>" + m + "</u></div>";
                for (const j of idx) {
                    const col = colordict[cts[j]] || '#000000';
                    const pct = (agg[j] * 100).toFixed(0);
                    h += "<div style='display:flex;align-items:center;margin:1px 0;'>"
                       + "<div style='width:32px;height:7px;background:#eee;margin-right:3px;border:1px solid #ccc;flex:none;'>"
                       + "<div style='width:" + Math.min(100, agg[j] * 100) + "%;height:7px;background:" + col + ";'></div></div>"
                       + "<span style='font-size:9px'>" + cts[j] + " " + pct + "%</span></div>";
                }
                return h;
            }
            function card(label, clusterList) {
                let h = "<div style='border:1px solid #bbb;border-radius:4px;padding:4px;width:165px;box-sizing:border-box;'>";
                h += "<div style='font-weight:bold;font-size:11px;border-bottom:1px solid #ccc;margin-bottom:2px;'>" + label + "</div>";
                for (const m of methods) { h += methodBlock(m, clusterList); }
                h += "</div>";
                return h;
            }

            // one card per selected cluster + an "Average" card, laid out as a grid
            let cards = "";
            for (const cl of active) {
                let spots = 0;
                for (const m of methods) { const c = comp['counts'][m][String(cl)]; if (c) { spots = c; break; } }
                cards += card("Cluster " + cl + " (" + spots + ")", [cl]);
            }
            if (active.length > 1) {
                cards += card("Average (" + active.join(', ') + ")", active);
            }
            div.text = "<div style='display:flex;flex-wrap:wrap;gap:6px;max-height:440px;overflow-y:auto;'>" + cards + "</div>";
        """)
        cluster_checkbox.js_on_change('active', comp_cb)

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
                    and the colors indicate the different clusters.
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
                In each spot, celltype proportions are represented with a pie chart
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
                    deviations among celltypes
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
            ">Spatial Transcriptomics Deconvolution Exploration  Interface</h1>
            <hr style="
                width: 80px;
                border: 2px solid #2c3e50;
                margin: 10px auto;
            ">
        </div>
        """
    # Create the Div widget
    info_box = Div(
        text= text1,
        width=700,
        height=200
    )

    # Update the button callbacks to refresh the Div text
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

    download_df = test_df.drop(columns=[col for col in test_df.columns if 'tooltip' in col])
    csv_source = ColumnDataSource({'data': [download_df.to_csv(index=False)]})
    download_button = Button(label="Download raw data", width=100, button_type='primary')
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

    p.add_tools(TapTool())

    leg = p.legend[0]
    leg.glyph_height = 20
    leg.glyph_width = 20
    leg.label_text_font_size = "15pt"
    leg.items = sorted(leg.items, key=lambda item: item.label['value'])
    p.add_layout(leg, 'right')

    title_box = Div(text=title_text, width=700, height=100)

    # Create the row for the buttons and slider
    buttons_col = column(
        show_all_button, *button_methods,
        rmsd_button, download_button , slider,
        cluster_controls,
        comp_div,
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
    p.legend.click_policy = "none"  # legend = color key; filtering is driven by the Cluster buttons
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

    # Average cell-type composition per cluster and per method (for the info panel)
    _clusters = pd.read_csv(data_clustered)
    _clusters = _clusters.set_index(_clusters.columns[0])['BayesSpace']
    comp_means, comp_counts, comp_celltypes = {}, {}, None
    for _m, _fp in zip(deconv_methods, norm_weights_filepaths):
        _props = pd.read_csv(_fp, sep='\t', index_col=0)
        _props.index.name = None
        if comp_celltypes is None:
            comp_celltypes = list(_props.columns)
        _df = _props.join(_clusters.rename('cluster'), how='inner').dropna(subset=['cluster'])
        _df['cluster'] = _df['cluster'].astype(int)
        _grp = _df.groupby('cluster')
        _means = _grp[comp_celltypes].mean()
        _cnt = _grp.size()
        comp_means[_m] = {str(int(c)): [float(v) for v in _means.loc[c].tolist()] for c in _means.index}
        comp_counts[_m] = {str(int(c)): int(_cnt.loc[c]) for c in _cnt.index}
    cluster_composition = {'celltypes': comp_celltypes, 'means': comp_means, 'counts': comp_counts}

    vis_with_separate_clusters_view(reduced_df=processed_data, image_path=image_path, deconv_methods=deconv_methods, nb_spots_samples=nb_spots_samples, n_largest_cell_types=n_largest_cell_types, output=output_html, cluster_composition=cluster_composition)
