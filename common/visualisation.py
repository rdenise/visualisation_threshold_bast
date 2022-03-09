import plotly.express as px
import plotly
import pandas as pd
import tqdm
import os

##########################################################################

def scatter2D_plotly(all_fam_file, name_tmp="tmp_interactive_scatter.html"):
    '''
    Function that take files that contains 5 columns:
    ['protein1', 'protein2', 'pident', 'evalue', 'coverage', 'fam']

    The function will plot the scatter plot in a temporary html file
    '''

    print('\n#####################')
    print('# Ploting the 2D plot')
    print('#####################\n')

    # Keep in memory the order of what we plot
    list_trace = []

    # Keep in memory the order of the seeds
    all_seeds = []

    # Help to keep an eye on the order of the dataframes plotted
    index = 0

    df_fam_all = pd.read_table(all_fam_file, usecols=['protein1', 'protein2','pident', 'evalue', 'coverage', 'fam'],
                            dtype = {'protein1': 'string',
                            'protein2': 'string',
                            'pident': 'float',
                            'evalue': 'float',
                            'coverage': 'float',
                            'fam': 'category',}
                                )
        
    df_fam_all['fam_name'] = df_fam_all.fam.apply(lambda x: x.split('family_')[-1])

    for fam, df_fam in tqdm.tqdm(df_fam_all.groupby('fam_name')) :

        # First element of family columns = in|out_family_seed
        seed = fam
        all_seeds.append(seed)
        
        # It is multiply by 6 because there is 6 traces for one plot 3 for in_fam et 3 for out_fam
        list_trace += [seed] * df_fam.fam.unique().shape[0] * 3
        
        # Change the name inside the columns to be more readable in the legend of the figure
        df_fam.replace(f"in_family_{seed}", 'Both in the family', inplace=True)
        df_fam.replace(f"out_family_{seed}", 'Only one in the family', inplace=True)
        
        # Reductio of the dataframe to remove the point in the same place on the plot 
        df_fam_drop = df_fam.drop_duplicates(['pident', 'coverage', 'fam']).reset_index(drop=True)  


        # create a figur efor the two histograms
        tmp_fig = px.scatter(df_fam, x='pident', y='coverage', color='fam', 
                             marginal_x='histogram', 
                             marginal_y='histogram',
                            color_discrete_sequence=px.colors.qualitative.Set1,
                            labels={"fam": "Pair of proteins"},
                            category_orders={"fam": ["Only one in the family", "Both in the family"]},
                            custom_data=['protein1','protein2', 'evalue'])

        # Create a figure for the scatter plot 
        tmp_fig_drop = px.scatter(df_fam_drop, x='pident', y='coverage', color='fam', 
                                 marginal_x='histogram', 
                                 marginal_y='histogram',
                                color_discrete_sequence=px.colors.qualitative.Set1,
                                labels={"fam": "Pair of proteins"},
                                category_orders={"fam": ["Only one in the family", "Both in the family"]},
                                custom_data=['protein1','protein2', 'evalue'],
                                )  

        # Update the information show when cliking on the point
        i = 0

        # To be more precise in the change we will change it manually:
        # trace 0 and 3 are the scatter plot
        # traco 1 and 4 are the histogram of the percentage of identity
        # trace 2 and 5 are the histogram of the coverage
        for data in tmp_fig.data:
            if i == 0 or i == 3 :
                # Replace the scatter plot of the drop dataframe instead of the full one
                data['customdata'] = tmp_fig_drop.data[i]['customdata']
                data['x'] = tmp_fig_drop.data[i]['x']
                data['y'] = tmp_fig_drop.data[i]['y']

                data['hovertemplate'] = "<br>".join([
                        "Protein 1 id: %{customdata[0]}",
                        "Protein 2 id: %{customdata[1]}",
                        "Percentage of identity: %{x}",
                        "Coverage: %{y}",
                        "E-value: %{customdata[2]}",
                    ])
            # Here only the value x (percentage id) and y (count)
            # are to change
            elif i == 1 or i == 4 :
                data['hovertemplate'] = "<br>".join([
                        "Percentage of identity: %{x}",
                        "Number: %{y}",
                    ])  
             # Here only the value y (coverage) and x (count)
            # are to change
            else :
                data['hovertemplate'] = "<br>".join([
                        "Coverage: %{y}",
                        "Number: %{x}",
                    ])              
            i+=1        

        for data in tmp_fig.data:
            # We want here to see the first figure but not the other one
            if index == 0 :
                data['visible'] = True
                fig = tmp_fig
            else :
                data['visible'] = False
                fig.add_trace(data)

        index += 1
        
    # Deal with the button to see one plot by seed 
    list_button = []

    # We create a list that will contains all the informations for the button
    # Here it is a button that when press will update the visualisation of the plot
    # So each buttons created will make visible the plot of the wanted seed
    for seed in all_seeds:
        list_button.append(dict(label=seed,
                method="update",
                args=[{"visible": [i == seed for i in list_trace]},
                    {"title": f"{seed} familly".capitalize(),
                     "title_x": 0.5,
                      },
                     ]
               ))

    # We update the layout of the figure to add the button created before
    fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=list(list_button),
                bordercolor = "#BEC8D9",  # default
                type="dropdown",  # "dropdown", "buttons"
                direction="down",
                showactive=True, # Highlights active dropdown item or active button if True.    
                y=0.9,
                x=-0.05
            )
        ])   


    # Set axis name
    fig.update_layout(
                    xaxis_title="Percentage of identity",
                    yaxis_title= "Coverage")

    # Set title
    fig.update_layout(
                     title={
                    'text': f"{all_seeds[0]} familly".capitalize(),
                    'font':{'size':30}
                     },)


    # Put the title in the middle
    fig.update_layout(title_x=0.5)

    plotly.offline.plot(fig, auto_open=False, filename=name_tmp)     

    return name_tmp

##########################################################################


def scatter3D_plotly(all_fam_file, name_tmp="tmp_interactive_scatter3D.html"):
    '''
    Function that take files that contains 5 columns:
    ['protein1', 'protein2', 'pident', 'evalue', 'coverage', 'fam']

    The function will plot the scatter plot in a temporary html file
    '''

    print('\n#####################')
    print('# Ploting the 3D plot')
    print('#####################\n')

    # Keep in memory the order of what we plot
    list_trace = []

    # Keep in memory the order of the seeds
    all_seeds = []

    # Help to keep an eye on the order of the dataframes plotted
    index = 0

    df_fam_all = pd.read_table(all_fam_file, usecols=['protein1', 'protein2','pident', 'evalue', 'coverage', 'fam'],
                            dtype = {'protein1': 'string',
                            'protein2': 'string',
                            'pident': 'float',
                            'evalue': 'float',
                            'coverage': 'float',
                            'fam': 'category',}
                                )

    df_fam_all['fam_name'] = df_fam_all.fam.apply(lambda x: x.split('family_')[-1])

    for fam, df_fam in tqdm.tqdm(df_fam_all.groupby('fam_name')) :
        
        df_fam = df_fam.drop_duplicates(['pident', 'coverage', 'evalue', 'fam']).reset_index(drop=True)

        # First element of family columns = in|out_family_seed
        seed = fam
        all_seeds.append(seed)

        # Change the name inside the columns to be more readable in the legend of the figure
        df_fam.replace(f"in_family_{seed}", 'Both in the family', inplace=True)
        df_fam.replace(f"out_family_{seed}", 'Only one in the family', inplace=True)    
        
        # list of the trace, two traces per seed (one in and one out)
        list_trace += [seed]*df_fam.fam.unique().shape[0]

        tmp_fig = px.scatter_3d(df_fam, x='pident', y='coverage', z='evalue',color='fam', 
                                color_discrete_sequence=px.colors.qualitative.Set1,
                               labels={"fam": "Pair of proteins"},
                               category_orders={"fam": ["Only one in the family", "Both in the family"]},
                               custom_data=['protein1','protein2', 'evalue'])
        
        tmp_fig.update_scenes(zaxis={'exponentformat':'e'})

        # Update the information show when cliking on the point
        tmp_fig.update_traces(
            hovertemplate="<br>".join([
                "Protein 1 id: %{customdata[0]}",
                "Protein 2 id: %{customdata[1]}",
                "Percentage of identity: %{x}",
                "Coverage: %{y}",
                "E-value: %{z}",
            ])
        )   
        
        for data in tmp_fig.data:
            if index == 0 :
                data['visible'] = False
                fig = tmp_fig
            else :
                data['visible'] = False
                fig.add_trace(data)
                
        index += 1
        
    # Deal with the button to see one plot by seed 
    list_button = []

    for seed in ['(no seed)'] + all_seeds:
        # We want here to see no figure but warn user of the size of it
        if seed == '(no seed)':
            list_button.append(dict(label=seed,
                    method="update",
                    args=[{"visible": [i == seed for i in list_trace]},
                          {"title": f"Choose a family to see the distribution.<br>WARNING:: Depending on your dataset you may need lot of memory",
                           "title_x": 0.5,
                          },
                         ]
                   ))
        else :
            list_button.append(dict(label=seed,
            method="update",
            args=[{"visible": [i == seed for i in list_trace]},
                  {"title": f"{seed} familly".capitalize(),
                           "title_x": 0.5,
                          },
                 ]
           ))

    # Deal with the scale of the zaxis from linear to log
    list_button2 = []

    for ztype in ['linear', 'log']:
            list_button2.append(dict(label=ztype,
            method="relayout",
            args=[{"scene.zaxis.type":ztype}]
           ))    


    # Update layout for the seeds
    fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=list(list_button),
                bordercolor = "#BEC8D9",  # default
                type="dropdown",  # "dropdown", "buttons"
                direction="down",
                showactive=True, # Highlights active dropdown item or active button if True.  
                y=0.8,
                x=-0.05
            ),
            dict(
                active=0,
                buttons=list(list_button2),
                bordercolor = "#BEC8D9",
                type="buttons",
                direction="left",
                showactive=True,
                y=0.9,
                x=-0.05
                ),
        ])   

    # Set axis name
    fig.update_layout(scene = dict(
                        xaxis_title='Percentage of identity',
                        yaxis_title='Coverage',
                        zaxis_title='E-value')
                     )
    # Set title
    fig.update_layout(
                     title={
                        'text': f"Choose a family of seed to see the distribution.<br>WARNING:: Depending on your dataset you may need lot of memory",
                        'font':{'size':30}
                     },
                     # annotate the buttons
                     annotations=[
                        dict(text="Evalue axis scale", x=-0.135, y=0.93,
                                             align="left", showarrow=False),
                        dict(text="Seeds", x=-0.13, y=0.83,
                                             showarrow=False),
                     ])

    # Put the title in the middle
    fig.update_layout(title_x=0.5)

    plotly.offline.plot(fig, auto_open=False, filename=name_tmp)    

    return name_tmp

##########################################################################


def fig2html(plot2D_file, plot3D_file, report, css):
    '''
    Function thst take the two html figure, read it and put it inside the 
    big report
    '''

    plot2D = ''
    plot3D = ''

    begin = False

    # Read plot line by line to remove non useful part of the html
    with open(plot2D_file, 'r', encoding='utf8') as r_file:
        for line in r_file:
            split_line = line.split()
            if line.strip().startswith('<script') :
                plot2D += line
                begin = True
            elif split_line[-1] == '</div>' and split_line[-2] == '</script>' :
                plot2D += line.replace('</div>', '')
                begin = False
            elif begin:
                plot2D += line
       
    # Read plot line by line to remove non useful part of the html 
    with open(plot3D_file, 'r', encoding='utf8') as r_file:
        for line in r_file :
            split_line = line.split()
            if line.strip().startswith('<script')  :
                plot3D += line
                begin = True
            elif split_line[-1] == '</div>' and split_line[-2] == '</script>' :
                plot3D += line.replace('</div>', '')
                begin = False
            elif begin:
                plot3D += line

    # Read the custom css to inject it directly to be portable
    with open(css, 'r', encoding='utf8') as r_file:
        css_string = r_file.read()

    # Write the html and inject css, plot2D and plot3D
    html_start = '''
    <!DOCTYPE html>
    <html lang="en">
        
        <head>
            <title>Distribition thresholds dashboard</title>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <meta http-equiv="X-UA-Compatible" content="IE=edge">
            
            <!-- Bootstrap CSS CDN -->
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">

            <!-- Our Custom CSS -->
            <style type="text/css" media="screen">
                ''' + css_string + '''
            </style>

            <!-- Font Awesome JS -->
            <script defer src="https://use.fontawesome.com/releases/v5.0.13/js/solid.js" integrity="sha384-tzzSw1/Vo+0N5UhStP3bvwWPq+uvzCMfrN1fEFe+xBmv1C/AtVX5K0uZtmcHitFZ" crossorigin="anonymous"></script>
            <script defer src="https://use.fontawesome.com/releases/v5.0.13/js/fontawesome.js" integrity="sha384-6OIrr52G08NpOFSZdxxz1xdNSndlD4vdcf/q2myIUVO0VsqaGHJsB0RaBE01VTOY" crossorigin="anonymous"></script>

        </head>
        
        <body>
        <div class="wrapper">
            <!-- Sidebar Holder -->
            <nav id="sidebar">
                <div class="sidebar-header">
                    <h3>Table of contents</h3>
                </div>

                <ul class="list-unstyled components">
                    <li>
                        <a href="#Plot2D">Plot2D: Percentage identity vs Coverage</a>
                    </li>
                    <li>
                        <a href="#Plot3D">Plot3D: Percentage identity vs Coverage vs E-value</a>
                    </li>
                </ul>
            </nav>

            <!-- Page Content Holder -->
            <div id="content">

                <nav class="navbar navbar-expand-lg navbar-light bg-light">
                    <div class="container-fluid">

                        <button type="button" id="sidebarCollapse" class="navbar-btn">
                            <span></span>
                            <span></span>
                            <span></span>
                        </button>
                        <button class="btn btn-dark d-inline-block d-lg-none ml-auto" type="button" data-toggle="collapse" data-target="#navbarSupportedContent" aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
                            <i class="fas fa-align-justify"></i>
                        </button>

                        <div class="collapse navbar-collapse" id="navbarSupportedContent">
                            <ul class="nav navbar-nav ml-auto">
                                <li class="nav-item">
                                    <a class="nav-link" href="#">Report</a>
                                </li>
                            </ul>
                        </div>
                    </div>
                </nav>
                <!-- Plot part : section name -->
                <h4><small>Plot 2D</small></h4>
                    <hr>
                    <!-- Title of the plot -->
                    <h2 id="Plot2D">Distribution of the family&apos;s hits : 2D plot</h2>
                        <!-- Description of the plot -->
                        <span class="fa fa-list-alt"></span> Description of the plot:
                         <ul>
                          <li>Family: This term is use to refer to your cluster and all the sequences belonging to this cluster using your thresholds</li>
                          <li>In family: Refer to pairs of proteins that are both inside the family</li>
                          <li>Out family: Refer to pairs of proteins with only one is inside the family</li>
                          <br>
                          This plot shows the population of pair of hits from the NCBI Blast ALL VS ALL part of the analysis. Each point represent a pair of 
                          two proteins. The hits are plot using the percentage of identity against the coverage of the blast alignment of the two proteins.
                          The histogram on top and bottom represent the distribution of respectively the percentage of identity and the coverage of your hits population.
                        </ul> 
                        <br>
                        <!-- Could be delete if your plot is not blue or red OR you can figure out how to change the color of the badges -->
                        <h5><span class="badge badge-danger">Pair containing protein outside of the family</span> <span class="badge badge-primary">Both proteins inside the family</span></h5>
                        <br>
                        <br>
                        <!-- In the CSS file there is a description of the class that allows to handle the height of the plot -->
                        <div class=plot_html>
                            <script type="text/javascript">window.PlotlyConfig = {MathJaxConfig: 'local'};</script>
                            ''' + plot2D + '''
                        </div>                
                        <br>
                        <br>
                <!-- Plot part : section name -->        
                <h4><small>Plot 3D</small></h4>
                    <hr>
                    <!-- Title of the plot -->
                    <h2 id="Plot3D">Distribution of the family&apos;s hits : 3D plot</h2>
                        <!-- Description of the plot -->
                        <span class="fa fa-list-alt"></span> Description of the plot:
                         <ul>
                          <li>Family: This term is use to refer to your cluster and all the sequences belonging to this cluster using your thresholds</li>
                          <li>In family: Refer to pairs of proteins that are both inside the family</li>
                          <li>Out family: Refer to pairs of proteins with only one is inside the family</li>
                          <br>
                          This plot shows the population of pair of hits from the NCBI Blast ALL VS ALL part of the analysis. Each point represent a pair of 
                          two proteins. The hits are plot using the percentage of identity against the coverage against the e-value of the blast alignment of 
                          the two proteins. You can change the scale of the e-value axis from linear to log by clicking on the button.
                          <br> 
                          <br>
                          <i class="fas fa-exclamation-triangle"></i> WARNING:: This plot is memory consuming be sure to have enough memory to handle it
                        </ul> 
                        <br>                    
                        <!-- Could be delete if your plot is not blue or red OR you can figure out how to change the color of the badges -->
                        <h5><span class="badge badge-danger">Pair containing protein outside of the family</span> <span class="badge badge-primary">Both proteins inside the family</span></h5>
                        <br>
                        <br>
                        <!-- In the CSS file there is a description of the class that allows to handle the height of the plot -->
                        <div class=plot_html>
                            <script type="text/javascript">window.PlotlyConfig = {MathJaxConfig: 'local'};</script>
                            ''' + plot3D + '''
                        </div>
                        <br>
                        <br>
                        <hr>
            </div>
            
            <!-- jQuery CDN - Slim version (=without AJAX) -->
            <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
            <!-- Popper.JS -->
            <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.0/umd/popper.min.js" integrity="sha384-cs/chFZiN24E4KMATLdqdvsezGxaGsi4hLGOzlXwp5UZB1LY//20VyM2taTB4QvJ" crossorigin="anonymous"></script>
            <!-- Bootstrap JS -->
            <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/js/bootstrap.min.js" integrity="sha384-uefMccjFJAIv6A+rW+L4AHf99KvxDjWSu1z9VI8SKNVmz4sk7buKt/6v9KI65qnm" crossorigin="anonymous"></script>

            <script type="text/javascript">
                $(document).ready(function () {
                    $('#sidebarCollapse').on('click', function () {
                        $('#sidebar').toggleClass('active');
                        $(this).toggleClass('active');
                    });
                });
            </script>
        </body>

    </html>


    '''

    # Write the report 
    with open(report, 'w') as w_file:
        w_file.write(html_start)

    os.remove(plot2D_file)
    os.remove(plot3D_file)

    return

##########################################################################