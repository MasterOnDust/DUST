import matplotlib.pyplot as plt
from . import maps 
import cartopy.crs as ccrs



def latex_plot():
    """
    create nice latex formating of plot. 
    """

    plt.rcParams['figure.figsize'] = (16, 16)
    plt.rcParams['font.size'] = 45
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['axes.linewidth'] = 1.5

    # --- Use latex-style on text in plots
    plt.rcParams['text.usetex'] = True

    # --- Custumize the length of the labels
    plt.rcParams["legend.labelspacing"] = 0.2
    plt.rcParams["legend.handlelength"] = 1.0
    plt.rcParams["legend.borderaxespad"] = 0.01

    # --- Ignore warnings for generated plot
    plt.rcParams.update({'figure.max_open_warning': 0})

    plt.linewidth=17.0
    font_size_plot = 20

    return font_size_plot

def base_figure_plot(nrows=1, ncols=1,terrain=True, **subplots_kw):
    """
    Setup the base matplotlib object for creating tracing the winds
    plot.  

    return matplotlib figure and cartopy.geoaxes

    """

    fig,ax =  plt.subplots(ncols=ncols,nrows=nrows,
        subplot_kw={'projection':ccrs.AlbersEqualArea(central_longitude=95.0, central_latitude=35.0)}, **subplots_kw)
    
    if terrain:
        maps.map_terrain_china()
    else:
        maps.map_china()
    ax.set_extent([75,115, 30, 48])

    return fig, ax