import cartopy as cr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.gridliner import Gridliner
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.colors as mcolors

from shapely.geometry import LineString, MultiLineString

from IPython import embed

def tracing_the_winds_map(ax):
    """
    DESCRIPTION:
    ============

        Generates a cartopy map of the Chinese Loess Plateu and the surrounding
        desert marked in, mountain ranges marked.

    USAGE:
    ======
        ax : cartopy.GeoAxes

        ax = tracing_the_winds_map(ax)

            returns: cartopy.geoaxes


    """

    shpfilename = shpreader.natural_earth(resolution='10m',
                                        category='physical',
                                        name='geography_regions_polys')

    boundary_10m = cfeature.NaturalEarthFeature('cultural',
                                            name ='admin_0_boundary_lines_land',
                                            scale ='10m',
                                            facecolor='none')
    prov_10 = shpreader.natural_earth(category = 'cultural',
                                            name ='admin_1_states_provinces',
                                            resolution= '10m')
    read_prov = shpreader.Reader(prov_10)

    prov_feature = read_prov.records()

    reader = shpreader.Reader(shpfilename)
    regionfeature = reader.records()
    China_desert = ['GOBI DESERT', 'Mu Us Desert', 'TAKLIMAKAN DESERT']
    China_mountain = ['TIAN SHAN', 'ALTUN MTS.', 'HIMALAYAS', 'PLATEAU OF TIBET']
    China_plateu = ['Loess Plateau']
    for feature in regionfeature: 
        if feature.attributes['name'] in China_desert:
            geom = [feature.geometry]

            # g = ax.add_geometries(geom, ccrs.PlateCarree(), facecolor='none', edgecolor='tan', linewidth=2)
            x = feature.geometry.centroid.x
            y = feature.geometry.centroid.y
            ax.text(x, y, feature.attributes['name'], color='darkorange', size=12, ha='center',
                    va='center', transform=ccrs.PlateCarree())
        elif feature.attributes['name'] in China_mountain:
            geom = [feature.geometry]

            g = ax.add_geometries(geom, ccrs.PlateCarree(), facecolor='none', edgecolor='slateblue', linewidth=1,
                                                        linestyle=':')
            x = feature.geometry.centroid.x
            y = feature.geometry.centroid.y
            ax.text(x, y, feature.attributes['name_en'], color='blue', size=12,
                    ha='center', va='center', transform=ccrs.PlateCarree())
        elif feature.attributes['name'] in China_plateu:
            geom = [feature.geometry]

            # g = ax.add_geometries(geom, ccrs.PlateCarree(), facecolor='none', edgecolor='peru', linewidth=2)
            x = feature.geometry.centroid.x
            y = feature.geometry.centroid.y
            ax.text(x, y, feature.attributes['name_en'], color='sienna', size=12,
                    ha='center', va='center', transform=ccrs.PlateCarree())

    for prov in prov_feature:
        if prov.attributes['admin'] == 'China' and prov.attributes['type_en'] == 'Province':

            ax.add_geometries([prov.geometry],crs =ccrs.PlateCarree(), facecolor='none', edgecolor='gray',
                            alpha=.8, linestyle=':')


    ax.add_feature(boundary_10m, edgecolor='gray')
    ax.coastlines('10m', color='gray', alpha=0.8)
    ax.set_extent([70,120, 25, 50], crs=ccrs.PlateCarree())

    gl = ax.gridlines(transform = ccrs.PlateCarree(), draw_labels = True, linestyle ='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

def map_terrain_china(ax):
    """
    DESCRIPTION
    ===========
    
        Draws a Stamen terrain map. Takes a cartopy.geoaxes as argument

    USAGE
    =====

        ax = map_terrain_china(ax)

        ax: cartopy.GeoAxes

        retrun ax cartopy.GeoAxes

    """
    stamen_terrain = cimgt.Stamen('terrain-background')


    ax = map_china(ax)
    ax.add_image(stamen_terrain, 7)
    ax.add_feature(cr.feature.RIVERS)
    ax.add_feature(cr.feature.LAKES)
    return ax



def map_china(ax, lakes_and_rivers=False):
    """
    DESCRIPTION
    ===========
    
        Draws a maps of china including administrative districts. 
        Takes a cartopy.geoaxes as argument

    USAGE
    =====

        ax = map_terrain_china(ax)

        ax: cartopy.GeoAxes

        retrun ax cartopy.GeoAxes

    """
    boundary_10m = cfeature.NaturalEarthFeature('cultural',
                                            name ='admin_0_boundary_lines_land',
                                            scale ='10m',
                                            facecolor='none')
    prov_10 = shpreader.natural_earth(category = 'cultural',
                                            name ='admin_1_states_provinces',
                                            resolution= '10m')

    read_prov = shpreader.Reader(prov_10)

    prov_feature = read_prov.records()

    for prov in prov_feature:
        if prov.attributes['admin'] == 'China' and prov.attributes['type_en'] == 'Province':
            ax.add_geometries([prov.geometry],crs =ccrs.PlateCarree(), facecolor='none', edgecolor='gray',
                            alpha=.8, linestyle=':')

    ax.add_feature(boundary_10m, edgecolor='gray')
    ax.coastlines('10m', color='gray', alpha=0.8)
    ax.set_extent([70,120, 25, 50], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linestyle ='--', xformatter=LONGITUDE_FORMATTER,
                   yformatter=LATITUDE_FORMATTER)
    gl.top_labels = False
    gl.right_labels = False
    if lakes_and_rivers:
        ax.add_feature(cr.feature.RIVERS)
        ax.add_feature(cr.feature.LAKES)
    return ax


def base_map_func(ax):
    """
    DESCRIPTION
    ===========
    
        Draw a basic map with coastlines and contry borders
        Takes cartop.GeoAxes instance as argument
    
    USAGE
    =====

        ax = base_map_func(ax)

        ax : cartopy.GeoAxes
    """
    
    ax.coastlines()
    ax.add_feature(cr.feature.BORDERS)
    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, color = 'grey', alpha = 0.6, linestyle = '--')
    gl.top_labels = False
    gl.right_labels = False

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax
