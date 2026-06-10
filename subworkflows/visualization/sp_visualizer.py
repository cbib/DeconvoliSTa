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

def process_data(norm_weights_filepath, st_coords_filepath, data_clustered, deconv_method, n_largest_cell_types, scale_factor, data_clustered2=None):
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

    # Keep the spot barcode as an explicit column so it survives the inter-method
    # merges in post_process_data (those reset the index); it is needed downstream to
    # align the full per-spot proportions used by the cell-type heatmap.
    merged_df["barcode"] = merged_df.index

    data_with_clusters = pd.read_csv(data_clustered)
    clusters_col =  pd.DataFrame(data_with_clusters["BayesSpace"]).set_index(data_with_clusters["Unnamed: 0"])
    merged_df["Cluster"] = clusters_col

    # Optional second clustering (for the Seurat <-> BayesSpace toggle)
    if data_clustered2 is not None:
        dwc2 = pd.read_csv(data_clustered2)
        merged_df["Cluster2"] = pd.DataFrame(dwc2["BayesSpace"]).set_index(dwc2["Unnamed: 0"])


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
    columns_of_interest = ['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres','Cluster' , "in_tissue"] + [f"{deconv_method}_Deconv_cell{i + 1}_norm_value" for i in range(n_largest_cell_types)] \
        + [f"{deconv_method}_Deconv_cell{i + 1}" for i in range(n_largest_cell_types)]
    if data_clustered2 is not None:
        columns_of_interest = columns_of_interest + ['Cluster2']
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


def post_process_data(norm_weights_filepaths, st_coords_filepath, data_clustered, deconv_methods, n_largest_cell_types, scale_factor, data_clustered2=None):
    norm_weights_dfs = [process_data(props, st_coords_filepath,data_clustered, deconv_methods[index], n_largest_cell_types, scale_factor = scale_factor, data_clustered2=data_clustered2)\
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
from bokeh.models import BoxAnnotation, Label, Plot, Rect, Text, Button, CustomJS, Div,Slider, PanTool, CheckboxButtonGroup, RadioButtonGroup
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from bokeh.layouts import column, row, gridplot,Spacer
from PIL import Image
import numpy as np

import math

    
   
def vis_with_separate_clusters_view(reduced_df, image_path, deconv_methods, nb_spots_samples, n_largest_cell_types, output, cluster_composition=None, clustering_labels=None, full_props=None, show_legend=False, show_figure=False):
    from bokeh.models import LinearColorMapper, ColorBar
    from bokeh.palettes import Viridis256

    image_display_infos = get_image_display_infos(image_path)
    image_display_infos = {x: int(np.ceil(image_display_infos[x]/2)) if x != "image_path" else image_display_infos[x] for x in image_display_infos}

    test_df = reduced_df[reduced_df["in_tissue"] == 1].head(nb_spots_samples).copy().reset_index(drop=True)
    n_spots = len(test_df)

    # Spot barcodes in plotting order, and the full (all cell types) per-spot proportion
    # matrix per method, aligned to that order. These feed the clickable cell-type heatmap:
    # clicking a type in the side panel recolors the active method's map by that type's
    # abundance. Canonical cell-type list = columns of the first method's proportions.
    barcodes = test_df['barcode'].tolist() if 'barcode' in test_df.columns else list(test_df.index)
    all_cell_types = list(full_props[deconv_methods[0]].columns) if full_props else []
    full_props_ordered = {}
    if full_props:
        for _m in deconv_methods:
            full_props_ordered[_m] = full_props[_m].reindex(index=barcodes, columns=all_cell_types)

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

    # --- Cluster plot: a single scatter colored by the ACTIVE clustering (toggleable). ---
    # Filtering is done via per-spot alpha (like the deconv plots), so the clustering toggle
    # only needs to swap the per-spot color/cluster fields.
    p = figure(width=900, height=700, title="Clustering results",
               x_axis_label='x', y_axis_label='y')
    p.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)

    has_two_clusterings = 'Cluster2' in test_df.columns
    clu0 = [int(c) for c in test_df['Cluster'].tolist()]
    clu1 = [int(c) for c in (test_df['Cluster2'] if has_two_clusterings else test_df['Cluster']).tolist()]

    # Color the ALTERNATIVE clustering so each of its domains keeps the color of the PRIMARY
    # domain it spatially overlaps most -> a region keeps its color when toggling Seurat<->BayesSpace.
    # The alternative usually has more clusters; the surplus domains (no primary left to match)
    # get fresh colors from the rest of the palette (not used by the primary).
    def _matched_colormap(primary_per_spot, alt_per_spot):
        from collections import Counter
        overlap = Counter(zip(alt_per_spot, primary_per_spot))   # (alt, primary) -> shared spots
        prim_ids = sorted(set(primary_per_spot))
        alt_ids = sorted(set(alt_per_spot))
        # greedy one-to-one assignment, strongest overlaps first
        alt_to_prim, used_prim = {}, set()
        for (a, pr), _cnt in sorted(overlap.items(), key=lambda kv: -kv[1]):
            if a in alt_to_prim or pr in used_prim:
                continue
            alt_to_prim[a] = pr
            used_prim.add(pr)
        used_colors = {clusters_colordict.get(p) for p in prim_ids}
        reserve = [c for _k, c in sorted(clusters_colordict.items()) if c not in used_colors]
        cmap, ri = {}, 0
        for a in alt_ids:
            if a in alt_to_prim:
                cmap[a] = clusters_colordict.get(alt_to_prim[a], '#000000')
            else:
                cmap[a] = reserve[ri] if ri < len(reserve) else '#000000'
                ri += 1
        return cmap

    primary_cmap = {c: clusters_colordict.get(c, '#000000') for c in set(clu0)}
    alt_cmap = _matched_colormap(clu0, clu1) if has_two_clusterings else primary_cmap
    col0 = [primary_cmap.get(c, '#000000') for c in clu0]
    col1 = [alt_cmap.get(c, '#000000') for c in clu1]
    cluster_source = ColumnDataSource(dict(
        x=all_x, y=all_y,
        cluster=clu0[:], cluster_0=clu0, cluster_1=clu1,
        color=col0[:], color_0=col0, color_1=col1,
        alpha=[1.0] * n_spots,
        tooltip_data=test_df['tooltip_data'].tolist()
    ))
    p.scatter(x='x', y='y', size=5, marker="circle", fill_color='color',
              fill_alpha='alpha', line_width=0, source=cluster_source)
    p.add_tools(HoverTool(tooltips="<div style='width:220px'>@tooltip_data{safe}</div>"))

    # cluster ids per clustering (index 0 = primary, 1 = alternative); active = primary
    cluster_ids_list = [sorted(set(clu0)), sorted(set(clu1))]
    cluster_ids = cluster_ids_list[0]

    # --- Deconv plots: vectorized — 1 source per method, n_largest_cell_types renderers ---
    # Replaces the previous n_spots × n_cell_types individual sources/renderers loop.
    deconv_plots = []
    deconv_sources = []
    # Per-method renderers for the cell-type heatmap: the pie wedges (default view) and a
    # single-circle layer colored by the selected cell type(s) abundance (heatmap view).
    wedge_renderers = []
    hm_renderers = []
    method_pies = []          # per method: K wedge renderers for the "pie of selected types" layer
    K_PIE = 8                 # max number of cell types shown in a per-spot pie
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
            'cluster': clu0[:],
            'cluster_0': clu0,
            'cluster_1': clu1,
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
        # Full per-spot proportions for every cell type (for the clickable heatmap), plus the
        # running heatmap value (sum of the currently selected cell types), updated in JS.
        if full_props:
            fpm = full_props_ordered[method]
            for ct in all_cell_types:
                source_data['prop_' + ct] = [
                    round(float(v), 4) if pd.notna(v) else 0.0 for v in fpm[ct].tolist()
                ]
            # Single cell-type map: per-spot color (the type's own color) + opacity (its
            # proportion) + hover text, all filled in by JS when a type is clicked.
            source_data['hm_color'] = ['#000000'] * n_spots
            source_data['hm_alpha'] = [0.0] * n_spots
            source_data['hm_tip'] = [''] * n_spots
            # "Pie of selected types" layer: K wedges per spot (angles renormalized among the
            # selected types), one shared opacity per spot. All filled in by JS.
            source_data['pie_alpha'] = [0.0] * n_spots
            source_data['pie_tip'] = [''] * n_spots
            for k in range(K_PIE):
                source_data[f'pie_start_{k}'] = [0.0] * n_spots
                source_data[f'pie_end_{k}'] = [0.0] * n_spots
                source_data[f'pie_color_{k}'] = ['#000000'] * n_spots
        shared_source = ColumnDataSource(source_data)

        plot = figure(width=900, height=700, title=f"Deconvolution results - {method}",
                      x_axis_label='x', y_axis_label='y',
                      x_range=p.x_range, y_range=p.y_range)
        plot.image_url(url='url', x='x', y='y', w='w', h='h', alpha='alpha', source=image_source)

        method_wedges = []
        for j in range(n_largest_cell_types):
            wr = plot.wedge(x='x', y='y', radius=4.7,
                            start_angle=f'start_{j}', end_angle=f'end_{j}',
                            fill_color=f'color_{j}', fill_alpha='alpha', line_width=0, source=shared_source)
            method_wedges.append(wr)

        # Single cell-type layer: each spot drawn in the cell type's own color, opacity ∝ its
        # proportion at that spot (a "where is this cell type" map). Hidden until a type is clicked.
        hm_renderer = plot.scatter(x='x', y='y', size=6, marker="circle",
                                   fill_color='hm_color', fill_alpha='hm_alpha',
                                   line_width=0, source=shared_source)
        hm_renderer.visible = False

        # Pie-of-selected-types layer: K wedge renderers, hidden until "pie mode" is on.
        pie_rends = []
        for k in range(K_PIE):
            pr = plot.wedge(x='x', y='y', radius=4.7,
                            start_angle=f'pie_start_{k}', end_angle=f'pie_end_{k}',
                            fill_color=f'pie_color_{k}', fill_alpha='pie_alpha',
                            line_width=0, source=shared_source)
            pr.visible = False
            pie_rends.append(pr)

        plot.add_tools(HoverTool(tooltips="<div style='width:220px'>@tooltip_data{safe}</div>",
                                 renderers=method_wedges))
        plot.add_tools(HoverTool(tooltips="<div style='width:220px'>@hm_tip{safe}</div>",
                                 renderers=[hm_renderer]))
        plot.add_tools(HoverTool(tooltips="<div style='width:240px'>@pie_tip{safe}</div>",
                                 renderers=pie_rends))
        plot.visible = False
        deconv_plots.append(plot)
        deconv_sources.append(shared_source)
        wedge_renderers.append(method_wedges)
        hm_renderers.append(hm_renderer)
        method_pies.append(pie_rends)

    # --- RMSD plot: single vectorized scatter (replaces n_spots individual sources) ---
    error_values = test_df["error_value"].tolist()
    min_val, max_val = min(error_values), max(error_values)
    color_map = LinearColorMapper(palette=Viridis256, low=min_val, high=max_val)

    rmsd_source = ColumnDataSource({
        'x': all_x, 'y': all_y,
        'cluster': clu0[:], 'cluster_0': clu0, 'cluster_1': clu1,
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
    # A cluster selection applies everywhere (like the zoom), via per-spot alpha. Select/Deselect all.
    # state_src holds: the active clustering index (0 = primary, 1 = alternative); the active
    # view ('clusters' | a method index as string | 'compare'); and the list of cell types
    # selected for the heatmap (stored as a single-cell array to keep CDS columns equal-length).
    state_src = ColumnDataSource(dict(clustering=[0], view=['clusters'], selected=[[]],
                                      expanded=[0], pie_mode=[0]))

    cluster_filter_code = """
        const idx = state_src.data['clustering'][0];
        const ids = cluster_ids_list[idx];
        const active = new Set(checkbox.active.map(i => ids[i]));
        const all_sources = deconv_sources.concat([rmsd_source, cluster_source]);
        for (const s of all_sources) {
            const cl = s.data['cluster'];
            const al = s.data['alpha'];
            for (let i = 0; i < cl.length; i++) { al[i] = active.has(cl[i]) ? 1 : 0; }
            s.change.emit();
        }
    """
    cluster_checkbox = CheckboxButtonGroup(
        labels=[str(c) for c in cluster_ids],
        active=list(range(len(cluster_ids))),
        width=300
    )
    cluster_checkbox.js_on_change('active', CustomJS(
        args=dict(checkbox=cluster_checkbox, deconv_sources=deconv_sources,
                  rmsd_source=rmsd_source, cluster_source=cluster_source,
                  cluster_ids_list=cluster_ids_list, state_src=state_src),
        code=cluster_filter_code))
    select_all_btn = Button(label="Select all clusters", width=120, button_type='success')
    select_all_btn.js_on_click(CustomJS(
        args=dict(checkbox=cluster_checkbox),
        code="checkbox.active = Array.from(Array(checkbox.labels.length).keys());"))
    deselect_all_btn = Button(label="Deselect all clusters", width=120, button_type='warning')
    deselect_all_btn.js_on_click(CustomJS(
        args=dict(checkbox=cluster_checkbox),
        code="checkbox.active = [];"))

    # Color key (legend) for the active clustering's clusters
    _card_style = ("background:#fffdf8;border:1px solid #e6e1d4;border-radius:10px;"
                   "padding:11px 13px;box-shadow:0 1px 3px rgba(20,48,43,.06);"
                   "font-family:'IBM Plex Sans',sans-serif;")
    _ttl_style = "font-weight:600;color:#14302b;font-size:13px;display:block;margin-bottom:8px;"

    def _legend_html(ids, cmap):
        items = "".join(
            "<div style='display:flex;align-items:center;margin:3px 0;'>"
            "<div style='width:13px;height:13px;border-radius:3px;background:" + cmap.get(c, '#000000') +
            ";border:1px solid rgba(0,0,0,.15);margin-right:8px;'></div>"
            "<span style='font-size:12px;color:#41463f'>Cluster " + str(c) + "</span></div>"
            for c in ids)
        return ("<div style='" + _card_style + "'><span style='" + _ttl_style + "'>Cluster colors</span>"
                + items + "</div>")
    legend_html = [_legend_html(cluster_ids_list[0], primary_cmap), _legend_html(cluster_ids_list[1], alt_cmap)]
    # Context-sensitive panel to the right of the map: the cluster color key on the Clusters
    # view, and the active method's clickable cell-type composition on a deconvolution view.
    # Default view is "clusters", so it starts as the cluster legend.
    side_div = Div(text=legend_html[0], width=248)

    cluster_controls_children = [
        Div(text="<b>Clusters</b> (filter applied to all views):", width=320),
        row(select_all_btn, deselect_all_btn),
        cluster_checkbox
    ]

    # Clustering toggle (only when a second clustering is provided): swaps color + cluster
    # fields everywhere, rebuilds the cluster buttons, and re-applies the filter/composition.
    if has_two_clusterings:
        clustering_toggle = RadioButtonGroup(labels=clustering_labels, active=0)
        clustering_toggle.js_on_change('active', CustomJS(
            args=dict(state_src=state_src, cluster_source=cluster_source,
                      deconv_sources=deconv_sources, rmsd_source=rmsd_source,
                      cluster_ids_list=cluster_ids_list, checkbox=cluster_checkbox),
            code="""
            const idx = cb_obj.active;
            state_src.data['clustering'] = [idx];
            const suf = '_' + idx;
            cluster_source.data['color'] = cluster_source.data['color' + suf].slice();
            cluster_source.data['cluster'] = cluster_source.data['cluster' + suf].slice();
            for (const s of deconv_sources.concat([rmsd_source])) {
                s.data['cluster'] = s.data['cluster' + suf].slice();
            }
            cluster_source.change.emit();
            const ids = cluster_ids_list[idx];
            checkbox.labels = ids.map(c => String(c));
            // resetting the checkbox triggers the cluster filter + side-panel re-render,
            // which refreshes the legend / composition for the new clustering.
            checkbox.active = Array.from(Array(ids.length).keys());
            """))
        cluster_controls_children = [Div(text="<b>Clustering:</b>", width=320),
                                     clustering_toggle] + cluster_controls_children

    cluster_controls = column(*cluster_controls_children)

    # --- Average cell-type composition panel (cards spread across the full width). ---
    # Initial content is rendered server-side (Python) so the tables show on load, without
    # needing a first interaction; the CustomJS below updates them on cluster/clustering changes.
    def _comp_card_html(comp, ids):
        cts = comp['celltypes']
        def method_block(m, cl_list):
            means = comp['means'].get(m, {}); counts = comp['counts'].get(m, {})
            agg = [0.0] * len(cts); tw = 0
            for cl in cl_list:
                key = str(cl)
                if key not in means: continue
                w = counts[key]; v = means[key]
                for j in range(len(cts)): agg[j] += w * v[j]
                tw += w
            if tw == 0: return ""
            agg = [a / tw for a in agg]
            order = sorted(range(len(cts)), key=lambda j: agg[j], reverse=True)[:4]
            h = "<div style='margin-top:4px;font-size:10px;color:#0f766e;font-weight:600'>" + m + "</div>"
            for j in order:
                col = colordict.get(cts[j], '#000000'); pct = "%.0f" % (agg[j] * 100)
                h += ("<div style='display:flex;align-items:center;margin:1.5px 0;'>"
                      "<div style='width:34px;height:7px;border-radius:3px;background:#efece2;margin-right:4px;flex:none;overflow:hidden;'>"
                      "<div style='width:" + str(min(100, agg[j] * 100)) + "%;height:7px;background:" + col + ";'></div></div>"
                      "<span style='font-size:9.5px;color:#41463f'>" + cts[j] + " <span class='pct'>" + pct + "%</span></span></div>")
            return h
        def card(label, cl_list):
            h = "<div class='dv-card' style='padding:7px 9px;width:168px;box-sizing:border-box;'>"
            h += "<div class='hd' style='font-size:11.5px;border-bottom:1px solid #ece7da;padding-bottom:3px;margin-bottom:2px;'>" + label + "</div>"
            for m in deconv_methods: h += method_block(m, cl_list)
            return h + "</div>"
        cards = ""
        for cl in ids:
            spots = 0
            for m in deconv_methods:
                c = comp['counts'].get(m, {}).get(str(cl))
                if c: spots = c; break
            cards += card("Cluster " + str(cl) + " (" + str(spots) + ")", [cl])
        if len(ids) > 1:
            cards += card("Average (" + ", ".join(str(c) for c in ids) + ")", ids)
        return "<div style='display:flex;flex-wrap:wrap;gap:8px;'>" + cards + "</div>"

    _comp_initial = (_comp_card_html(cluster_composition[0], cluster_ids_list[0])
                     if cluster_composition is not None
                     else "<i>Select one or more cluster(s) to see their cell-type composition.</i>")
    comp_div = Div(text=_comp_initial, sizing_mode='stretch_width', height=360)
    if cluster_composition is not None:
        comp_cb = CustomJS(args=dict(checkbox=cluster_checkbox, cluster_ids_list=cluster_ids_list,
                                     comp_list=cluster_composition, div=comp_div,
                                     methods=deconv_methods, colordict=colordict, state_src=state_src),
            code="""
            const idx = state_src.data['clustering'][0];
            const ids = cluster_ids_list[idx];
            const comp = comp_list[idx];
            const active = checkbox.active.map(i => ids[i]);
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
                let h = "<div style='margin-top:4px;font-size:10px;color:#0f766e;font-weight:600'>" + m + "</div>";
                for (const j of idx) {
                    const col = colordict[cts[j]] || '#000000';
                    const pct = (agg[j] * 100).toFixed(0);
                    h += "<div style='display:flex;align-items:center;margin:1.5px 0;'>"
                       + "<div style='width:34px;height:7px;border-radius:3px;background:#efece2;margin-right:4px;flex:none;overflow:hidden;'>"
                       + "<div style='width:" + Math.min(100, agg[j] * 100) + "%;height:7px;background:" + col + ";'></div></div>"
                       + "<span style='font-size:9.5px;color:#41463f'>" + cts[j] + " <span class='pct'>" + pct + "%</span></span></div>";
                }
                return h;
            }
            function card(label, clusterList) {
                let h = "<div class='dv-card' style='padding:7px 9px;width:168px;box-sizing:border-box;'>";
                h += "<div class='hd' style='font-size:11.5px;border-bottom:1px solid #ece7da;padding-bottom:3px;margin-bottom:2px;'>" + label + "</div>";
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
            div.text = "<div style='display:flex;flex-wrap:wrap;gap:8px;'>" + cards + "</div>";
        """)
        cluster_checkbox.js_on_change('active', comp_cb)

    # --- Right-of-map interactive panel (legend on Clusters view; clickable composition on a
    # method view). The rendering logic + cell-type toggle handlers live on `window` so the
    # HTML rows (built inside a Bokeh Div) can call back into the Bokeh models. They are defined
    # once (guarded by window.dvInit), capturing the model references from the first callback
    # that runs; every view/filter callback re-runs window.dvRender(). ---
    _compare_note = ("<div style='" + _card_style + "'><span style='" + _ttl_style + "'>Method comparison</span>"
                     "<div style='font-size:11.5px;color:#41463f;line-height:1.5'>Spots are colored by "
                     "how much the methods disagree (Viridis scale, right of the map). Brighter = less "
                     "agreement between methods.</div></div>")
    side_init_js = """
    if (!window.dvInit) {
        window.dvState = state_src;
        window.dvM = {side_div, legend_html, comp_list, cluster_ids_list, checkbox, methods,
                      colordict, deconv_sources, hm_renderers, wedge_renderers, compare_note,
                      clustering_labels, method_pies, K_PIE};
        window.dvToggleCT = function(name) {
            const sel = window.dvState.data['selected'][0];
            const k = sel.indexOf(name);
            if (k >= 0) sel.splice(k, 1); else sel.push(name);
            window.dvRender();
        };
        window.dvClearCT = function() { window.dvState.data['selected'][0] = []; window.dvRender(); };
        window.dvToggleOthers = function() {
            const e = window.dvState.data['expanded']; e[0] = e[0] ? 0 : 1; window.dvRender();
        };
        window.dvTogglePie = function() {
            const e = window.dvState.data['pie_mode']; e[0] = e[0] ? 0 : 1; window.dvRender();
        };
        window.dvRender = function() {
            const M = window.dvM, S = window.dvState.data;
            const cidx = S['clustering'][0], view = S['view'][0], selected = S['selected'][0];
            const ids = M.cluster_ids_list[cidx];
            const activeClusters = M.checkbox.active.map(i => ids[i]);
            // reset every method to its default (per-spot deconvolution pies) view
            for (let mi = 0; mi < M.methods.length; mi++) {
                M.hm_renderers[mi].visible = false;
                for (const w of M.wedge_renderers[mi]) w.visible = true;
                for (const pr of M.method_pies[mi]) pr.visible = false;
            }
            if (view === 'clusters') { M.side_div.text = M.legend_html[cidx]; return; }
            if (view === 'compare')  { M.side_div.text = M.compare_note; return; }

            const mi = parseInt(view), m = M.methods[mi];
            const comp = M.comp_list[cidx], cts = comp['celltypes'];
            const means = comp['means'][m], counts = comp['counts'][m];
            let agg = new Array(cts.length).fill(0), tw = 0;
            for (const cl of activeClusters) {
                const key = String(cl);
                if (!(key in means)) continue;
                const w = counts[key], v = means[key];
                for (let j = 0; j < cts.length; j++) agg[j] += w * v[j];
                tw += w;
            }
            if (tw > 0) for (let j = 0; j < cts.length; j++) agg[j] /= tw;
            const order = Array.from(cts.keys()).sort((a, b) => agg[b] - agg[a]);
            const selSet = new Set(selected);
            const expanded = S['expanded'][0], TOPN = 12;

            // All styling is INLINE (Bokeh renders each Div in its own shadow root, so global
            // CSS classes would not reach the panel — that is why the swatches went missing).
            const cardCss = "background:#fffdf8;border:1px solid #e6e1d4;border-radius:10px;"
                + "padding:11px 13px;box-shadow:0 1px 3px rgba(20,48,43,.06);font-family:'IBM Plex Sans',sans-serif;";
            const clab = (M.clustering_labels && M.clustering_labels[cidx]) ? M.clustering_labels[cidx] : '';
            let html = "<div style='" + cardCss + "'>"
                + "<div style='font-weight:600;color:#14302b;font-size:13px;'>" + m + " &mdash; composition</div>"
                + "<div style='font-size:10.5px;color:#6b7069;margin:2px 0 8px'>" + clab + " &middot; "
                + activeClusters.length + "/" + ids.length + " domains &middot; click a type to see its spots</div>";
            if (selected.length)
                html += "<div onclick='window.dvClearCT()' style='display:inline-block;font-size:11px;"
                      + "font-weight:600;color:#b45309;cursor:pointer;margin-bottom:8px'>&times; clear "
                      + selected.length + " selected</div>";
            // With >=2 types selected, offer a per-spot pie of the selected types (vs the
            // dominant-type map). The pie shows their relative balance; opacity = their total.
            if (selected.length >= 2) {
                const pm = S['pie_mode'][0];
                html += "<div onclick='window.dvTogglePie()' style='display:inline-block;font-size:11px;"
                      + "font-weight:600;color:#0f766e;cursor:pointer;margin:0 0 8px 12px'>"
                      + (pm ? "&#9632; back to dominant-type map" : "&#9685; create pie chart from selected")
                      + "</div>";
            }

            const rowCss = "display:flex;align-items:center;gap:8px;padding:4px 6px;border-radius:6px;cursor:pointer;";
            function ctRow(j) {
                const name = cts[j], col = M.colordict[name] || '#000000';
                const pct = (agg[j] * 100).toFixed(1), on = selSet.has(name);
                const hl = on ? "background:#eaf2f0;box-shadow:inset 0 0 0 1.5px #0f766e;" : "";
                return "<div onclick=\\"window.dvToggleCT('" + name + "')\\" style='" + rowCss + hl + "'>"
                     + "<span style='width:13px;height:13px;border-radius:2px;border:1px solid rgba(0,0,0,.2);"
                     + "flex:none;background:" + col + "'></span>"
                     + "<span style='flex:1;min-width:0;font-size:11.5px;color:#2f342e;white-space:nowrap;"
                     + "overflow:hidden;text-overflow:ellipsis'>" + name + "</span>"
                     + "<span style='font-family:monospace;font-size:11px;color:#0f766e;flex:none'>" + pct + "%</span></div>";
            }
            // Always show the selected types, plus the top TOPN; the rest hide behind "others".
            const head = order.slice(0, TOPN);
            const tail = order.slice(TOPN);
            let rows = "";
            for (const j of head) rows += ctRow(j);
            // selected types that fell into the tail stay visible even when collapsed
            if (!expanded) for (const j of tail) { if (selSet.has(cts[j])) rows += ctRow(j); }
            if (expanded) for (const j of tail) rows += ctRow(j);
            html += "<div style='max-height:540px;overflow-y:auto;margin-top:2px'>" + rows + "</div>";
            if (tail.length > 0) {
                const lbl = expanded ? "&#9650; show top " + TOPN + " only"
                                     : "&#9660; show " + tail.length + " more cell types";
                html += "<div onclick='window.dvToggleOthers()' style='font-size:11px;color:#0f766e;"
                      + "cursor:pointer;font-weight:600;margin-top:7px;padding:3px 2px'>" + lbl + "</div>";
            }
            M.side_div.text = html + "</div>";

            const pieMode = S['pie_mode'][0] && selected.length >= 2;
            if (pieMode) {
                // Per-spot pie of the selected types: wedge angles renormalized among the
                // selected types (relative balance), one opacity per spot = their total
                // proportion (scaled to the strongest visible spot) -> faint where they are scarce.
                const src = M.deconv_sources[mi], d = src.data, al = d['alpha'];
                const n = al.length, TAU = 2 * Math.PI;
                const nsel = Math.min(selected.length, M.K_PIE), sel = selected.slice(0, nsel);
                for (let k = 0; k < M.K_PIE; k++) {
                    const cArr = d['pie_color_' + k], c = (k < nsel) ? (M.colordict[sel[k]] || '#000000') : '#000000';
                    for (let i = 0; i < n; i++) cArr[i] = c;
                }
                const sums = new Array(n); let hi = 0;
                for (let i = 0; i < n; i++) {
                    let s = 0; const vals = new Array(nsel);
                    for (let k = 0; k < nsel; k++) { const arr = d['prop_' + sel[k]]; const v = arr ? arr[i] : 0; vals[k] = v; s += v; }
                    sums[i] = s;
                    if (al[i] > 0 && s > hi) hi = s;
                    let ang = 0;
                    for (let k = 0; k < nsel; k++) {
                        const st = d['pie_start_' + k], en = d['pie_end_' + k];
                        if (s > 0) { st[i] = ang * TAU; ang += vals[k] / s; en[i] = ang * TAU; } else { st[i] = 0; en[i] = 0; }
                    }
                    for (let k = nsel; k < M.K_PIE; k++) { d['pie_start_' + k][i] = 0; d['pie_end_' + k][i] = 0; }
                }
                if (hi <= 0) hi = 1;
                const pa = d['pie_alpha'], pt = d['pie_tip'];
                for (let i = 0; i < n; i++) {
                    pa[i] = al[i] > 0 ? Math.min(1, sums[i] / hi) : 0;
                    if (al[i] > 0 && sums[i] > 0) {
                        let tip = "<b>selected types</b> &mdash; " + (sums[i] * 100).toFixed(1) + "% of this spot<br>";
                        for (let k = 0; k < nsel; k++) {
                            const arr = d['prop_' + sel[k]], v = arr ? arr[i] : 0, c = M.colordict[sel[k]] || '#000000';
                            tip += "<div style='display:flex;align-items:center;'><div style='width:9px;height:9px;background:"
                                 + c + ";margin-right:5px;'></div>" + sel[k] + ": " + (v * 100).toFixed(1) + "%</div>";
                        }
                        pt[i] = tip;
                    } else pt[i] = '';
                }
                src.change.emit();
                M.hm_renderers[mi].visible = false;
                for (const w of M.wedge_renderers[mi]) w.visible = false;
                for (let k = 0; k < M.K_PIE; k++) M.method_pies[mi][k].visible = (k < nsel);
            } else if (selected.length > 0) {
                // Spot map: each spot shown in the dominant selected type's OWN color, with opacity
                // proportional to that type's proportion (scaled to the strongest visible spot).
                // -> "where is this cell type" rather than a rainbow heatmap.
                const src = M.deconv_sources[mi], d = src.data;
                const hc = d['hm_color'], ha = d['hm_alpha'], ht = d['hm_tip'], al = d['alpha'];
                const best = new Array(al.length), bestv = new Array(al.length);
                let hi = 0;
                for (let i = 0; i < al.length; i++) {
                    let bj = null, bv = 0;
                    for (const name of selected) { const arr = d['prop_' + name]; const v = arr ? arr[i] : 0; if (v > bv) { bv = v; bj = name; } }
                    best[i] = bj; bestv[i] = bv;
                    if (al[i] > 0 && bv > hi) hi = bv;
                }
                if (hi <= 0) hi = 1;
                for (let i = 0; i < al.length; i++) {
                    if (best[i] === null) { hc[i] = '#000000'; ha[i] = 0; ht[i] = ''; continue; }
                    const name = best[i], col = M.colordict[name] || '#000000';
                    hc[i] = col;
                    ha[i] = al[i] > 0 ? Math.min(1, bestv[i] / hi) : 0;
                    ht[i] = "<div style='display:flex;align-items:center;'>"
                          + "<div style='width:10px;height:10px;background:" + col + ";margin-right:5px;'></div>"
                          + "<span style='color:blue;'>" + name + "</span>: " + (bestv[i] * 100).toFixed(1) + "%</div>";
                }
                src.change.emit();
                M.hm_renderers[mi].visible = true;
                for (const w of M.wedge_renderers[mi]) w.visible = false;
            }
        };
        window.dvInit = true;
    }
    """
    _side_args = dict(state_src=state_src, side_div=side_div, legend_html=legend_html,
                      comp_list=(cluster_composition if cluster_composition is not None else []),
                      cluster_ids_list=cluster_ids_list, checkbox=cluster_checkbox,
                      methods=deconv_methods, colordict=colordict, deconv_sources=deconv_sources,
                      hm_renderers=hm_renderers, wedge_renderers=wedge_renderers,
                      compare_note=_compare_note, method_pies=method_pies, K_PIE=K_PIE,
                      clustering_labels=(clustering_labels if clustering_labels else ['Primary', 'Alternative']))

    def make_side_cb(set_view_js=""):
        return CustomJS(args=dict(_side_args), code=side_init_js + set_view_js + "\n window.dvRender();")

    # Re-render the side panel whenever the cluster filter changes (keeps the per-method
    # average + heatmap scaling in sync with the selected clusters).
    cluster_checkbox.js_on_change('active', make_side_cb())

    text1 = """<div class="dv-panel"><h3>Clusters view</h3>
        Each point is a spot, colored by its cluster (spatial domain). Use the
        <b>Clustering</b> toggle and the <b>Cluster</b> buttons to filter every view.</div>"""
    text2 = """<div class="dv-panel"><h3>Deconvolution by method</h3>
        Each spot is a pie chart of its cell-type proportions for the selected method.</div>"""
    text3 = """<div class="dv-panel"><h3>Method comparison</h3>
        Spots colored by disagreement between methods &mdash; RMSD for two methods, or the
        standard deviation of per-cell-type standard deviations for more than two.</div>"""
    title_text = """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,500;9..144,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@500&display=swap');
    html, body, .bk-root { background:#f3f1ea !important; }
    .bk-root { color:#1f2421; }
    .dv-header { padding:4px 2px 2px; }
    .dv-header .sub { font-family:'IBM Plex Mono',monospace; font-size:10.5px; letter-spacing:0.14em; text-transform:uppercase; color:#0f766e; }
    .dv-header h1 { font-family:'Fraunces',serif; font-weight:600; font-size:24px; line-height:1.1; letter-spacing:-0.01em; color:#14302b; margin:6px 0 0; }
    .dv-header .rule { height:3px; width:52px; background:#0f766e; margin-top:11px; border-radius:2px; }
    .dv-panel { font-family:'IBM Plex Sans',sans-serif; background:#fffdf8; border:1px solid #e6e1d4; border-left:3px solid #0f766e;
                border-radius:10px; padding:11px 14px; box-shadow:0 1px 3px rgba(20,48,43,.06); font-size:13px; line-height:1.55; color:#41463f; }
    .dv-panel h3 { font-family:'Fraunces',serif; font-size:14.5px; margin:0 0 5px; color:#14302b; font-weight:600; }
    .dv-legend { font-family:'IBM Plex Sans',sans-serif; background:#fffdf8; border:1px solid #e6e1d4; border-radius:10px;
                 padding:11px 13px; box-shadow:0 1px 3px rgba(20,48,43,.06); }
    .dv-legend .ttl { font-family:'Fraunces',serif; font-weight:600; color:#14302b; font-size:13px; display:block; margin-bottom:7px; }
    .dv-card { font-family:'IBM Plex Sans',sans-serif; background:#fffdf8; border:1px solid #e6e1d4; border-radius:9px;
               box-shadow:0 1px 2px rgba(20,48,43,.06); }
    .dv-card .hd { font-family:'Fraunces',serif; font-weight:600; color:#14302b; }
    .dv-card .pct { font-family:'IBM Plex Mono',monospace; color:#0f766e; }
    .dv-clear { display:inline-block; font-size:11px; font-weight:600; color:#b45309; cursor:pointer; margin:0 0 8px; }
    .dv-clear:hover { text-decoration:underline; }
    .dv-ctlist { max-height:564px; overflow-y:auto; margin-right:-5px; padding-right:5px; }
    .dv-ctlist::-webkit-scrollbar { width:7px; }
    .dv-ctlist::-webkit-scrollbar-thumb { background:#d8d2c2; border-radius:4px; }
    .dv-ctrow { display:flex; align-items:center; gap:8px; padding:4px 6px; border-radius:6px; cursor:pointer; transition:background .12s; }
    .dv-ctrow:hover { background:#f0ede3; }
    .dv-ctrow.on { background:#eaf2f0; box-shadow:inset 0 0 0 1.5px #0f766e; }
    .dv-ctrow .sw { width:13px; height:13px; border-radius:2px; border:1px solid rgba(0,0,0,.2); flex:none; }
    .dv-ctrow .nm { font-size:11.5px; color:#2f342e; flex:1; min-width:0; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
    .dv-ctrow .pc { font-family:'IBM Plex Mono',monospace; font-size:11px; color:#0f766e; flex:none; }
    .bk-btn-primary  { background:#14302b !important; border-color:#14302b !important; color:#f3f1ea !important; }
    .bk-btn-primary:hover { background:#1f4a42 !important; }
    .bk-btn-success  { background:#0f766e !important; border-color:#0f766e !important; color:#f3f1ea !important; }
    .bk-btn-warning  { background:#b45309 !important; border-color:#b45309 !important; color:#f7efe2 !important; }
    .bk-btn { font-family:'IBM Plex Sans',sans-serif !important; border-radius:7px !important; font-weight:500 !important; }
    </style>
    <div class="dv-header">
        <div class="sub">CBiB &middot; spatial transcriptomics</div>
        <h1>Deconvolution Exploration Interface</h1>
        <div class="rule"></div>
    </div>
    """
    # Create the Div widget
    info_box = Div(
        text= text1,
        width=360,
        height=130
    )

    # Update the button callbacks to refresh the Div text
    show_all_button = Button(label="Show Clusters", width=100, button_type = 'primary')
    show_all_button.js_on_click(CustomJS(args=dict(p=p,  rmsd_plot = rmsd_plot, deconv_plots=deconv_plots, info_box=info_box, text1 = text1), code="""
        p.visible = true;
        deconv_plots.forEach((p) => {p.visible = false});
        rmsd_plot.visible = false;
        info_box.text = text1;
    """))
    # Side panel: switch to the cluster legend on the Clusters view.
    show_all_button.js_on_click(make_side_cb("window.dvState.data['view'] = ['clusters'];"))
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
        # Side panel: show this method's clickable cell-type composition.
        button.js_on_click(make_side_cb(f"window.dvState.data['view'] = ['{index}'];"))
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
    # Side panel: comparison note (no cluster legend on the compare view).
    rmsd_button.js_on_click(make_side_cb("window.dvState.data['view'] = ['compare'];"))

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

    # (no per-cluster Bokeh legend: the cluster plot is a single scatter colored by the
    #  active clustering; the Cluster buttons act as the legend/control.)

    title_box = Div(text=title_text, width=1100, height=84)

    # Left column: controls only (view buttons, slider, clustering toggle, cluster filter)
    buttons_col = column(
        show_all_button, *button_methods,
        rmsd_button, download_button , slider,
        cluster_controls,
        Spacer(height=20)
    )

    # Sidebar: info + controls (no header, no composition)
    controls_column = column(
        info_box,
        buttons_col,
    )
    # Center: just the active plot
    center_column = column(
        p,
        *deconv_plots,
        rmsd_plot,
    )
    # Main row: [sidebar | plot | context panel]. The right panel shows the cluster color key
    # on the Clusters view, and the active method's clickable cell-type composition otherwise.
    main_row = row(
        controls_column,
        center_column,
        side_div,
    )
    # Full-width: header on top, composition as a bottom band (cards spread across the width)
    layout = column(
        title_box,
        main_row,
        Spacer(height=10),
        comp_div,
        sizing_mode='stretch_width'
    )

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
    # Optional: a second clustering file (argv[10]) + button labels (argv[11], e.g. "Seurat,BayesSpace")
    data_clustered2 = argv[10] if len(argv) > 10 and argv[10] not in ('', 'none', 'None') else None
    clustering_labels = argv[11].split(',') if len(argv) > 11 else ["Primary", "Alternative"]

    print("Processing data ...\n")
    processed_data = post_process_data(norm_weights_filepaths=norm_weights_filepaths, st_coords_filepath=st_coords_filepath, data_clustered=data_clustered, \
                                 deconv_methods=deconv_methods, n_largest_cell_types=n_largest_cell_types, scale_factor=scale_factor, data_clustered2=data_clustered2)

    print(f"Deconvolution methods {deconv_methods}")
    nb_spots_samples = processed_data.shape[0]
    print(f"Generating vis with {nb_spots_samples} spots and top {n_largest_cell_types} cells...\n")

    # Average cell-type composition per cluster and per method (for the info panel), per clustering.
    def _composition(clustering_file):
        clu = pd.read_csv(clustering_file)
        clu = clu.set_index(clu.columns[0])['BayesSpace']
        means, counts, celltypes = {}, {}, None
        for _m, _fp in zip(deconv_methods, norm_weights_filepaths):
            _props = pd.read_csv(_fp, sep='\t', index_col=0)
            _props.index.name = None
            if celltypes is None:
                celltypes = list(_props.columns)
            _df = _props.join(clu.rename('cluster'), how='inner').dropna(subset=['cluster'])
            _df['cluster'] = _df['cluster'].astype(int)
            _grp = _df.groupby('cluster')
            _means = _grp[celltypes].mean()
            _cnt = _grp.size()
            means[_m] = {str(int(c)): [float(v) for v in _means.loc[c].tolist()] for c in _means.index}
            counts[_m] = {str(int(c)): int(_cnt.loc[c]) for c in _cnt.index}
        return {'celltypes': celltypes, 'means': means, 'counts': counts}

    comp_primary = _composition(data_clustered)
    comp_alt = _composition(data_clustered2) if data_clustered2 is not None else comp_primary
    cluster_composition = [comp_primary, comp_alt]

    # Full (all cell types) per-spot proportions per method, for the clickable cell-type heatmap.
    full_props = {}
    for _m, _fp in zip(deconv_methods, norm_weights_filepaths):
        _fpdf = pd.read_csv(_fp, sep='\t', index_col=0)
        _fpdf.index.name = None
        full_props[_m] = _fpdf

    vis_with_separate_clusters_view(reduced_df=processed_data, image_path=image_path, deconv_methods=deconv_methods, nb_spots_samples=nb_spots_samples, n_largest_cell_types=n_largest_cell_types, output=output_html, cluster_composition=cluster_composition, clustering_labels=clustering_labels, full_props=full_props)
