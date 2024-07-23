#!/usr/bin/env python
# coding: utf-8

"""
Created on Mon Apr 22 16:46:25 2024

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

colordict = {
    "malignantcell": "#99CB66",
    "macrophage": "#FF3366",  
    "muralcell": "#996600",  
    "dendriticcell": "#00FFFF",  
    "microglialcell": "#CCCC99",  
    "monocyte": "#FF9999", 
    "oligodendrocyte": "#666666",
    "endothelialcell": "#00FFCC", 
    "matureTcell": "#0066FF",  
    "oligodendrocyteprecursorcell": "#CCFF99", 
    "mastcell": "#FF00CC",
    "Bcell": "#F0E442",
    "plasmacell": "#669900",
    "naturalkillercell": "#FFE6E6",  
    "astrocyte": "#00FFCC",
    "radialglialcell": "#009933",  
    "neuron": "#6600FF"
}


# <! ------------------------------------------------------------------------!>
# <!                           DATA PREPARATION                              !>
# <! ------------------------------------------------------------------------!>
from PIL import Image
def process_data(norm_weights_filepath, st_coords_filepath, data_clustered, image_path, n_largest_cell_types, scale_factor):
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
    im = Image.open(image_path).convert("RGB")
    # Merge coordinate df and cell weight df
    image_display_infos = {
        "image_path" : image_path,
        "x0" : 0,
        "y0" : 0,
        "im_w" : im.size[0],
        "im_h" : im.size[1],
        
    }
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
        merged_df[''.join(['Deconv_cell', str(i + 1)])] = cell_type_storage_arrays[i]
        merged_df[''.join(['Deconv_cell', str(i + 1), '_value'])] = cell_value_storage_arrays[i]

    # Since we only consider the top N cell types, we need to correct the weight
    # values so that the scatterpies account to the totality of the circle (sum of weights == 1)

    deconv_weight_columns = [f"Deconv_cell{i + 1}_value" for i in range(n_largest_cell_types)]

    # Create new normalized columns
    for i in range(n_largest_cell_types):

        # Calculate the sum of the top cell type weights
        total = merged_df.loc[:, deconv_weight_columns].sum(axis=1)

        # Create column with corrected weight values
        merged_df[''.join(['Deconv_cell', str(i + 1), '_norm_value'])] =  merged_df[''.join(['Deconv_cell', str(i + 1), '_value'])] / total


    # SLim down the df by selecting columns of interest only
    columns_of_interest = ['pxl_row_in_fullres', 'pxl_col_in_fullres','Cluster' , "in_tissue"] + [f"Deconv_cell{i + 1}_norm_value" for i in range(n_largest_cell_types)] \
        + [f"Deconv_cell{i + 1}" for i in range(n_largest_cell_types)]
    reduced_df = merged_df.loc[:, columns_of_interest]
    return reduced_df, image_display_infos


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
def vis_with_separate_clusters_view(reduced_df, image_display_infos, nb_spots_samples, output , show_legend = False, show_figure = False ):
        image_display_infos = {x: int(np.ceil(image_display_infos[x]/2)) if x != "image_path" else image_display_infos[x] for x in image_display_infos}
        # Smaller sample
        test_df = reduced_df[reduced_df["in_tissue"] == 1].head(nb_spots_samples).copy()
        # Create a single tooltip column for each circle
        test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'.join( \
                                                [f"<span style='color: red;'> Spot</span> : (x = { row['pxl_col_in_fullres']:.2f}, y = {-row['pxl_row_in_fullres']:.2f})"] ),\
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
            url=[ image_display_infos.get("image_path")],
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

        # Create a single tooltip column for each circle
        test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'.join([
            f"<div style='display:flex;align-items:center;'>"
            f"<div style='width:10px;height:10px;background-color:{colordict.get(row[f'Deconv_cell{i+1}'], '#000000')};margin-right:5px;'></div>"
            f"<span style='color: blue;'>{row[f'Deconv_cell{i+1}']}</span>: {row[f'Deconv_cell{i+1}_norm_value']*100:.2f}%"
            f"</div>"
            for i in range(n_largest_cell_types)
        ] +  [f"<span style='color: red;'> Spot</span> : (x = {row['pxl_col_in_fullres']:.2f}, y = {-row['pxl_row_in_fullres']:.2f})"]), axis=1)
        data["tooltip_data"] = test_df['tooltip_data'].tolist()
        for i in range(1, n_largest_cell_types + 1):
            data[f'DeconvCell{i}'] = test_df[f'Deconv_cell{i}'].tolist()
            data[f'DeconvCell{i}_w'] = test_df[f'Deconv_cell{i}_norm_value'].tolist()
        # Convert dictionary to dataframe
        df = pd.DataFrame(data)
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
                wedge = plot.wedge(x='x', y='y', radius=50,
                        start_angle=start_angle, end_angle=end_angle,
                        line_width=0, fill_color=colors[i],
                        legend_label=f"Cluster {row['Cluster']}", source=circle_source, visible=False)
                start_angle = end_angle
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
        show_all_button.js_on_click(CustomJS(args=dict(p=p, plot=plot, info_box=info_box, text1 = text1), code="""
            p.visible = true;
            plot.visible = false;
            info_box.text = text1;
        """))
        spacer = Spacer(width=50)  # Adjust the width as needed
        button = Button(label="Show deconvolution by Cluster", width=120)
        button.js_on_click(CustomJS(args=dict(p=p, plot=plot, info_box=info_box, text2 = text2), code="""
            p.visible = false;
            plot.visible = true;
            info_box.text = text2;
        """))
        spacer1 = Spacer(width=100)
        spacer2 = Spacer(width=100)
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
        tap_tool = TapTool()
        p.add_tools(tap_tool)
        # Associez le callback au TapTool
        hover = HoverTool(tooltips="""
            <div style="width:220px">
                @tooltip_data
            </div>
        """)
        p.add_tools(hover)
        plot.add_tools(hover)
        leg = p.legend[0]
        leg.glyph_height = 40
        leg.glyph_width = 40
        leg.label_text_font_size = "20pt"
        p.add_layout(leg,'left')
        leg_plot = plot.legend[0]
        leg_plot.glyph_width = 0
        plot.add_layout(leg_plot,'left')

        from bokeh.layouts import column, row # to avoid a conflict with row from pandas
        # Créez le layout
        buttons_row = row(show_all_button, spacer, button, spacer1, download_button, spacer2, slider)
        layout = column(buttons_row, info_box, p, plot)

        p.legend.location = "top_right"
        p.legend.click_policy = "hide"
        p.legend.visible = True
        plot.legend.location = "top_right"
        plot.legend.click_policy = "hide"
        plot.visible = False
        plot.legend.visible = True

        if show_figure:
            show(layout)
        output_file(output, mode='inline')
        save(layout)






# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>
def vis_with_proportions(reduced_df, nb_spots_samples, output, show_legend = False, width =1000, height = 1000 , show_fig = False):
        # Smaller sample
        test_df = reduced_df.iloc[1:nb_spots_samples, ].copy()
        # Create a single tooltip column for each circle
        test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'.join([
            f"<div style='display:flex;align-items:center;'>"
            f"<div style='width:10px;height:10px;background-color:{colordict.get(row[f'Deconv_cell{i+1}'], '#000000')};margin-right:5px;'></div>"
            f"<span style='color: blue;'>{row[f'Deconv_cell{i+1}']}</span>: {row[f'Deconv_cell{i+1}_norm_value']*100:.2f}%"
            f"</div>"
            for i in range(n_largest_cell_types)
        ] + [f"<span style='color: red;'> Spot</span>{ row['x'], row['y'] }"]), axis=1)

        # Update the data dictionary
        data = {
            'x': [y / 100 for y in test_df.y.tolist()],
            'y': [-x / 100 for x in test_df.x.tolist()],
            'x_full': test_df.y.tolist(),
            'y_full': [-x for x in test_df.x.tolist()],
            'tooltip_data': test_df['tooltip_data'].tolist(),
        }

        for i in range(1, n_largest_cell_types + 1):
            data[f'DeconvCell{i}'] = test_df[f'Deconv_cell{i}'].tolist()
            data[f'DeconvCell{i}_w'] = test_df[f'Deconv_cell{i}_norm_value'].tolist()
        # Convert dictionary to dataframe
        df = pd.DataFrame(data)
        # Convert dataframe to a ColumnDataSource
        # Initialize the Bokeh plot
        p = figure(width=width, height=height,
                   title="Deconvolution results",
                   x_axis_label='x',
                   y_axis_label='y',
                   output_backend="webgl",
                  )
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
                wedge = p.wedge(x='x', y='y', radius=0.05,
                        start_angle=start_angle, end_angle=end_angle,
                        line_color="white", fill_color=colors[i],
                        legend_label=f"{cell_types[i]}", source=circle_source)

                start_angle = end_angle
        leg = p.legend[0]
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

        if show_fig :
            show(layout)
        output_file(output)
        save(layout)
        return layout

# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>

def vis_type_majoritaires(reduced_df, nb_spots_samples, output, show_figure = False, show_legend = False, width =1000, height = 1000 ):
        # Smaller df for testing
        test_df = reduced_df.iloc[1:nb_spots_samples, ].copy()
        # Create a single tooltip column for each circle
        test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'. \
                                                      join([f"<span style='color: blue;'>{row[f'Deconv_cell{1}']}</span>: {row[f'Deconv_cell{1}_norm_value']*100:.2f}%"] \
                                                           + [f"<span style='color: red;'> Spot</span> {row['x']:.2f} , {row['y']:.2f}"] ),\
                                                      axis=1)
        # Update the data dictionary
        data = {
            'x': [y / 100 for y in test_df.y.tolist()],
            'y': [-x / 100 for x in test_df.x.tolist()],
            'x_full': test_df.y.tolist(),
            'y_full': [-x for x in test_df.x.tolist()],
            'tooltip_data': test_df['tooltip_data'].tolist(),
        }
        data[f'DeconvCell{1}'] = test_df[f'Deconv_cell{1}'].tolist()
        data[f'DeconvCell{1}_w'] = test_df[f'Deconv_cell{1}_norm_value'].tolist()
        # Convert dictionary to dataframe
        df = pd.DataFrame(data)
        # Initialize the Bokeh plot
        p = figure(width =1000, height = 1000,
                    title = "Deconvolution results",
                    x_axis_label = 'x',
                    y_axis_label = 'y',
                   output_backend="webgl"
                    )
        for index, row in df.iterrows():
            x, y = row['x'], row['y']
            categorie = row[f'DeconvCell{1}_w' ]
            cell_type = row[f'DeconvCell{1}']
            color = colordict[cell_type]
            # Create a single ColumnDataSource for all wedges in this circle
            circle_source = ColumnDataSource({
                'x': [x],
                'y': [y],
                'tooltip_data': [row['tooltip_data']]
            })
            scatter = p.scatter(x='x', y='y', size=15, 
                                marker="circle",  # Specify the marker shape
                                fill_color=color, line_color=color
                                , line_width=0, 
                                source=circle_source,
                                legend_label=f"{cell_type}")
        hover = HoverTool(tooltips="""
            <div style="width:200px">
                <h3>Proportions:</h3>
                @tooltip_data
            </div>
        """)
        # Add the hover tool to the plot
        p.add_tools(hover)
        l = p.legend[0]
        p.add_layout(l, 'left')
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
                    Cette vue montre le type cellulaire majoritaire de chaque spot"
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

        if show_figure:
            show(layout)
        output_file(output)
        save(layout)



# <! ------------------------------------------------------------------------!>
# <!                       BOKEH VISUALIZATION                               !>
# <! ------------------------------------------------------------------------!>


def vis_with_clustring(reduced_df, nb_spots_samples,  show_figure = False, show_legend = False, width =1000, height = 1000 ):
        # Smaller df for testing
        test_df = reduced_df.iloc[1:nb_spots_samples, ].copy()
        # Create a single tooltip column for each circle
        test_df['tooltip_data'] = test_df.apply(lambda row: '<br>'.join( \
                                                [f"<span style='color: red;'> Spot</span> : (x = { row['x']:.2f}, y = {row['y']:.2f})"] \
                                                + [f"<span style='color: blue;'> Cluster</span> { row['Cluster']}"] ),\
                                                axis=1)
        # Update the data dictionary
        data = {
            'x': [y / 100 for y in test_df.y.tolist()],
            'y': [-x / 100 for x in test_df.x.tolist()],
            'tooltip_data': test_df['tooltip_data'].tolist(),
            'Cluster' : test_df['Cluster'].tolist()
        }
        # Convert dictionary to dataframe
        df = pd.DataFrame(data)
        # Initialize the Bokeh plot
        p = figure(width =width, height = height,
                    title = "Deconvolution results",
                    x_axis_label = 'x',
                    y_axis_label = 'y',
                   output_backend="webgl"
                    )
        for index, row in df.iterrows():
            x, y = row['x'], row['y']
            cluster = row['Cluster']
            color = clusters_colordict[cluster]
            # Create a single ColumnDataSource for all wedges in this circle
            circle_source = ColumnDataSource({
                'x': [x],
                'y': [y],
                'tooltip_data': [row['tooltip_data']]
            })

            wedge = p.wedge(x='x', y='y', radius=0.02,
                    start_angle=0, end_angle=2*pi,
                    line_color="white", fill_color= color,  line_width = 0
                    , source=circle_source, legend_label= f"Cluster {cluster}")


        # Show no legend
        p.legend.visible= show_legend
        hover = HoverTool(tooltips="""
            <div style="width:200px">
                @tooltip_data
            </div>
        """)
        # Add the hover tool to the plot
        p.add_tools(hover)
        # Configurer la légende
        p.legend.location = "top_right"
        p.legend.click_policy = "hide"
        if show_figure:
            show(p)



if __name_ == "__main__":
    import os
    # norm_weights_filepath = "res_rctd_cluster/proportions_rctd_sample2"
    # st_coords_filepath = "tissue_positions_list_248.csv"
    # data_clustered = "seurat_metadata_UKF248_T_ST.csv"
    # n_largest_cell_types = 5
    # scale_factor = 0.24414062
    # nb_spots_samples = processed_data[0].shape[0] 
    # output_html = "visium_plot_sample2.html"

    argv = sys.argv
    norm_weights_filepath = argv[1]
    st_coords_filepath = argv[2]
    data_clustered = argv[3]
    n_largest_cell_types = argv[4]
    scale_factor = argv[5]
    nb_spots_samples = argv[6]
    output_html = argv[7]
    
    processed_data = process_data(norm_weights_filepath, st_coords_filepath,data_clustered, "original_tissue_images/tissue_hires_image_248.png", n_largest_cell_types, scale_factor = scale_factor)
    vis_with_separate_clusters_view(reduced_df=processed_data[0],image_display_infos= processed_data[1], nb_spots_samples = nb_spots_samples, output= output_html )

