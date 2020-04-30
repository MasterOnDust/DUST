import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import matplotlib.patheffects as PathEffects
from shapely.geometry import LineString, MultiLineString
import cartopy as cr
import matplotlib.colors as mcolors



"""
DESCRIPTION:
============

    Generates a cartopy map of the Chinese Loess Plateu and the surrounding 
    desert marked in, mountain ranges marked.  

USAGE:
======
    fig, ax = tracing_the_winds_map()

        returns:
             matplotlib.figure
             cartopy.geoaxes  
        

"""

def tracing_the_winds_map(figsize = (12,10)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection=ccrs.PlateCarree())

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
            
            g = ax.add_geometries(geom, ccrs.PlateCarree(), facecolor='none', edgecolor='tan', linewidth=2)
            x = feature.geometry.centroid.x        
            y = feature.geometry.centroid.y
            ax.text(x, y, feature.attributes['name_en'], color='darkorange', size=12, ha='center',
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
            
            g = ax.add_geometries(geom, ccrs.PlateCarree(), facecolor='none', edgecolor='peru', linewidth=2)
            x = feature.geometry.centroid.x        
            y = feature.geometry.centroid.y
            ax.text(x, y, feature.attributes['name_en'], color='sienna', size=12, 
                    ha='center', va='center', transform=ccrs.PlateCarree())
            
    for prov in prov_feature:
        if prov.attributes['admin'] == 'China' and prov.attributes['type_en'] == 'Province':
            

            
            ax.add_geometries(prov.geometry,crs =ccrs.PlateCarree(), facecolor='none', edgecolor='gray', 
                            alpha=.8, linestyle=':')

            
    ax.add_feature(boundary_10m, edgecolor='gray')
    ax.coastlines('10m', color='gray', alpha=0.8)
    ax.set_extent([70,120, 25, 50], crs=ccrs.PlateCarree())

    gl = ax.gridlines(transform = ccrs.PlateCarree(), draw_labels = True, linestyle ='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    return fig, ax

def base_map_func():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    # ax.add_feature(land_50m)
    ax.add_feature(cr.feature.BORDERS)
    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, color = 'grey', alpha = 0.6, linestyle = '--')
    gl.xlabels_top = False; gl.ylabels_right = False

    return ax