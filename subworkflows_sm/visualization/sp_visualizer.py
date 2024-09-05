# -*- coding: utf-8 -*-
"""
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
from bokeh.plotting import figure, output_file, save,output_notebook, show, curdoc

from bokeh.transform import cumsum
from bokeh.models import ColumnDataSource, HoverTool,Range1d
from bokeh.palettes import Category20
from bokeh.models import CustomJS, TapTool

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
from PIL import Image
def process_data(norm_weights_filepath, st_coords_filepath, data_clustered, deconv_method, n_largest_cell_types, scale_factor):
    # Read spatial deconvolution result CSV file
    norm_weights_df = pd.read_csv(norm_weights_filepath, sep = '\t')
    norm_weights_df.index.name = None
    # print(norm_weights_df.head())

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
    # print(max_weights)
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
def vis_with_separate_clusters_view(reduced_df, image_path, deconv_methods, nb_spots_samples, output , show_legend = False, show_figure = False ):
        image_display_infos = get_image_display_infos(image_path)
        image_display_infos = {x: int(np.ceil(image_display_infos[x]/2)) if x != "image_path" else image_display_infos[x] for x in image_display_infos}
        # Smaller sample
        test_df = reduced_df[reduced_df["in_tissue"] == 1].head(nb_spots_samples).copy()
        # Create a single tooltip column for each circle
        test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'.join( \
                                                [f"<span style='color: red;'> Spot</span> : (x = { row['pxl_col_in_fullres']/2:.2f}, y = {-row['pxl_row_in_fullres']/2:.2f})"] ),\
                                                axis=1)
        # Update the data dictionary
        data = {
            'x': [y/2 for y in test_df.pxl_col_in_fullres.tolist()],
            'y': [-x/2 for x in test_df.pxl_row_in_fullres.tolist()],
            'tooltip_data': test_df['tooltip_data'].tolist(),
            'Cluster' : test_df['Cluster'].tolist() ,
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

            scatter = p.scatter(x='x', y='y', size=7,
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
                  wedge = plot.wedge(x='x', y='y', radius=5,
                          start_angle=start_angle, end_angle=end_angle,
                          line_width=0, fill_color=colors[i],
                          legend_label=f"Cluster {row['Cluster']}", source=circle_source, visible=False)
                  start_angle = end_angle
          deconv_plots.append(plot)
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
                ">Information sur la visualisation</h3>
                <p style="margin-bottom: 0;">
                    Cette vue montre les clusters. Chaque point représente un spot,
                    et les couleurs indiquent les différents clusters.
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
                ">Information sur la visualisation</h3>
                <p style="margin-bottom: 0;">
                    Cette vue montre la déconvolution par cluster. Les graphiques représentent
                    la distribution des cellules dans chaque spot."
                </p>
            </div>
            """
        # Créez le widget Div
        info_box = Div(
            text= text1,
            width=1000,
            height=120
        )

        # Modifiez les callbacks des boutons pour mettre à jour le texte du Div
        show_all_button = Button(label="Show Clusters", width=100)
        show_all_button.js_on_click(CustomJS(args=dict(p=p, deconv_plots=deconv_plots, info_box=info_box, text1 = text1), code="""
            p.visible = true;
            deconv_plots.forEach((p) => {p.visible = false});
            info_box.text = text1;
        """))
        spacer = Spacer(width=50)  # Adjust the width as needed
        button_methods = []
        for index, method in enumerate(deconv_methods):
          button = Button(label=f"Show deconvolution by Cluster with {method}", width=120)
          button.js_on_click(CustomJS(args=dict(p=p, deconv_plots=deconv_plots, index = index, info_box=info_box, text2 = text2), code="""
              p.visible = false;
              deconv_plots.forEach((p) => {p.visible = false});
              deconv_plots[index].visible = true;
              info_box.text = text2;
          """))
          button_methods.append(button)
          spacer1 = Spacer(width=150)  # Adjust the width as needed
          button_methods.append(spacer1)

        spacer1 = Spacer(width=100)
        spacer2 = Spacer(width=100)
        # Assuming you have your data in a pandas DataFrame called 'df'
        csv_source = ColumnDataSource({'data': [df.drop(columns=[f"{method}_tooltip_data" for method in deconv_methods]).to_csv(index=False)]})
        download_button = Button(label="Download raw data", width=100)
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
          plot.add_layout(leg_plot,'left')
          plot.legend.location = "top_right"
          plot.legend.click_policy = "hide"
          plot.visible = False
          plot.legend.visible = True
        leg = p.legend[0]
        leg.glyph_height = 40
        leg.glyph_width = 40
        leg.label_text_font_size = "20pt"
        p.add_layout(leg,'left')



        from bokeh.layouts import column, row # to avoid a conflict with row from pandas
        # Créez le layout
        buttons_row = row(show_all_button, spacer, *button_methods, download_button, spacer2, slider)
        layout = column(buttons_row, info_box, p, *deconv_plots)

        p.legend.location = "top_right"
        p.legend.click_policy = "hide"
        p.legend.visible = True

        if show_figure:
            show(layout)
        output_file(output, mode='inline')
        save(layout)









# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>
# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>

def vis_with_proportions(reduced_df, image_display_infos, nb_spots_samples, output , show_legend = False, show_figure = False ):
        image_display_infos = {x: int(np.ceil(image_display_infos[x]/2)) if x != "image_path" else image_display_infos[x] for x in image_display_infos}
        # Smaller sample
        test_df = reduced_df[reduced_df["in_tissue"] == 1].head(nb_spots_samples).copy()
        # Create a single tooltip column for each circle
        test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'.join([
            f"<div style='display:flex;align-items:center;'>"
            f"<div style='width:10px;height:10px;background-color:{colordict.get(row[f'Deconv_cell{i+1}'], '#000000')};margin-right:5px;'></div>"
            f"<span style='color: blue;'>{row[f'Deconv_cell{i+1}']}</span>: {row[f'Deconv_cell{i+1}_norm_value']*100:.2f}%"
            f"</div>"
            for i in range(n_largest_cell_types)
        ] + [f"<span style='color: red;'> Spot</span> : (x = { row['pxl_col_in_fullres']:.2f}, y = {-row['pxl_row_in_fullres']:.2f})"] ), axis=1)

        # Update the data dictionary
        data = {
            'x': [y/2 for y in test_df.pxl_col_in_fullres.tolist()],
            'y': [-x/2 for x in test_df.pxl_row_in_fullres.tolist()],
            'tooltip_data': test_df['tooltip_data'].tolist(),
            'Cluster' : test_df['Cluster'].tolist() ,
        }
        for i in range(1, n_largest_cell_types + 1):
            data[f'DeconvCell{i}'] = test_df[f'Deconv_cell{i}'].tolist()
            data[f'DeconvCell{i}_w'] = test_df[f'Deconv_cell{i}_norm_value'].tolist()
        # Convert dictionary to dataframe
        df = pd.DataFrame(data)
        # Convert dataframe to a ColumnDataSource
        # Initialize the Bokeh plot
        p = figure(width = image_display_infos.get("im_w"), height = image_display_infos.get("im_h"),
                    title = "Deconvolution results",
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
    
        # Create a Div for displaying the message
        for index, row in df.iterrows():
            x, y = row['x'], row['y']
            categories = row[[f'DeconvCell{i+1}_w' for i in range(n_largest_cell_types)]].values
            cell_types = row[[f'DeconvCell{i+1}' for i in range(n_largest_cell_types)]].values
            colors = tuple([colordict[x] for x in cell_types])
            # Create a single ColumnDataSource for all wedges in this circle
            circle_source = ColumnDataSource({
                'x': [x],
                'y': [y],
                'tooltip_data': [row['tooltip_data']]
            })
            start_angle = 0
            for i, category_value in enumerate(categories):
                end_angle = start_angle + category_value * 2 * pi
                wedge = p.wedge(x='x', y='y', radius=5,
                        start_angle=start_angle, end_angle=end_angle,
                        line_width = 0, fill_color=colors[i],
                        legend_label=f"{cell_types[i]}", source=circle_source)

                start_angle = end_angle
        leg = p.legend[0]
        leg.glyph_height = 30
        leg.glyph_width = 30
        leg.label_text_font_size = "17pt"
        p.add_layout(leg, 'left')
        # # Show no legend
        p.legend.visible= show_legend
        hover = HoverTool(tooltips="""
            <div style="width:200px">
                <h3>Proportions:</h3>
                @tooltip_data
            </div>
        """)
        # Add the hover tool to the plot
        p.add_tools(hover)
        # Assuming you have your data in a pandas DataFrame called 'df'
        csv_source = ColumnDataSource({'data': [df.drop(columns=['tooltip_data']).to_csv(index=False)]})
        download_button = Button(label="Download raw data", width=100)
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
    
        text = """
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
                ">Information sur la visualisation</h3>
                <p style="margin-bottom: 0;">
                    Cette vue les resultats de la deconvolution, chaque element de la figure 
                    est un spot."
                </p>
            </div>
            """
        # Créez le widget Div
        info_box = Div(
            text= text,
            width=1000,
            height=120
        )
        from bokeh.layouts import column, row # to avoid a conflict with row from pandas

        layout = column(download_button, info_box, p)

        if show_figure :
            show(layout)
        output_file(output, mode = "inline")
        save(layout)
        return layout
# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>
if __name__ == "__main__":
    import sys
    argv = sys.argv
    norm_weights_filepaths = argv[1].split(',')
    st_coords_filepath = argv[2]
    data_clustered = argv[3]
    image_path = argv[4]
    n_largest_cell_types = int(argv[5])
    scale_factor = float(argv[6])
    output_html = argv[7]
    print("Processing data ...\n")
    # norm_weights_filepaths = ["drive/MyDrive/proportions_rctd_sample243_chunk1",  "drive/MyDrive/proportions_cell2location_UKF243_T_ST_1_raw_001_chunk_1.tsv"]
    norm_weights_dfs = [process_data(props, st_coords_filepath,data_clustered, deconv_methods[index], n_largest_cell_types, scale_factor = scale_factor)\
                        for index, props  in enumerate(norm_weights_filepaths)]
    processed_data = norm_weights_dfs[0]
    for i in range(1, len(norm_weights_dfs)):
        processed_data = pd.merge(processed_data, norm_weights_dfs[i])
    print(f"Deconvolution methods {deconv_methods}")
    deconv_methods = argv[8].split(',')
    print(f"Generating vis with {processed_data.shape[0]} spots and top {n_largest_cell_types} cells...\n")
    nb_spots_samples = processed_data.shape[0]
    vis_with_separate_clusters_view(reduced_df=processed_data, image_path = image_path, deconv_methods = deconv_methods,nb_spots_samples = nb_spots_samples, output= output_html )
