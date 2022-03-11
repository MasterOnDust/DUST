import numpy as np
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib.ticker import LogFormatterSciNotation


def _gen_log_clevs(dat_min, dat_max):
    """Creates a logarithmic color scale."""

    if dat_max > 0:
        dmx = int(np.round(np.log10(dat_max)))
    else:
        dmx = 1

    # TODO: What's the default value of dmn?
    if dat_min > 0:
        dmn = int(np.round(np.log10(dat_min)))
    elif dat_min == 0.0 or np.isnan(dat_min):
        dmn = dmx - 3

    # create equally spaced range
    # ERROR: dmn could be uninitialized
    if dmx == dmn:
        dmx = dmn + 1
    clevs = np.logspace(dmn, dmx, 100)

    return clevs


def _add_colorbar(im, cticks=None, label="", fmt="%d", log=False):
    """Adds Colorbar Nicely to figure"""
    ax = im.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(
        "right", size="3.5%", pad=0.05, axes_class=mpl.pyplot.Axes
    )
    if log:
        formatter = LogFormatterSciNotation(10, labelOnlyBase=False)
    else:
        formatter = fmt

    if isinstance(cticks, np.ndarray):
        fig.colorbar(im, ax=ax, cax=cax, label=label, format=formatter, ticks=cticks)
    else:
        fig.colorbar(im, ax=ax, cax=cax, label=label, format=formatter)


def _gen_flexpart_colormap(ctbfile=None, colors=None):
    """Generate the ast colormap for FLEXPART."""

    if ctbfile:
        try:
            colors = np.loadtxt(ctbfile)
        except:
            print("WARNING: cannot load ctbfile. using colors")
    if colors:
        name = "user_colormap"
    else:
        # AST Colorset for FLEXPART
        colors = [
            1.0000000e00,
            1.0000000e00,
            1.0000000e00,
            9.9607843e-01,
            9.1372549e-01,
            1.0000000e00,
            9.8431373e-01,
            8.2352941e-01,
            1.0000000e00,
            9.6470588e-01,
            7.1764706e-01,
            1.0000000e00,
            9.3333333e-01,
            6.0000000e-01,
            1.0000000e00,
            8.9019608e-01,
            4.4705882e-01,
            1.0000000e00,
            8.3137255e-01,
            2.0000000e-01,
            1.0000000e00,
            7.5686275e-01,
            0.0000000e00,
            1.0000000e00,
            6.6274510e-01,
            0.0000000e00,
            1.0000000e00,
            5.4901961e-01,
            0.0000000e00,
            1.0000000e00,
            4.0784314e-01,
            0.0000000e00,
            1.0000000e00,
            2.4705882e-01,
            0.0000000e00,
            1.0000000e00,
            7.4509804e-02,
            0.0000000e00,
            1.0000000e00,
            0.0000000e00,
            2.8235294e-01,
            1.0000000e00,
            0.0000000e00,
            4.8627451e-01,
            1.0000000e00,
            0.0000000e00,
            6.3137255e-01,
            1.0000000e00,
            0.0000000e00,
            7.4509804e-01,
            1.0000000e00,
            0.0000000e00,
            8.4705882e-01,
            1.0000000e00,
            0.0000000e00,
            9.3725490e-01,
            1.0000000e00,
            0.0000000e00,
            1.0000000e00,
            9.7647059e-01,
            0.0000000e00,
            1.0000000e00,
            8.9411765e-01,
            0.0000000e00,
            1.0000000e00,
            8.0000000e-01,
            0.0000000e00,
            1.0000000e00,
            6.9019608e-01,
            0.0000000e00,
            1.0000000e00,
            5.6470588e-01,
            0.0000000e00,
            1.0000000e00,
            4.0000000e-01,
            0.0000000e00,
            1.0000000e00,
            0.0000000e00,
            3.9607843e-01,
            1.0000000e00,
            0.0000000e00,
            5.6470588e-01,
            1.0000000e00,
            0.0000000e00,
            6.9019608e-01,
            1.0000000e00,
            0.0000000e00,
            7.9607843e-01,
            1.0000000e00,
            0.0000000e00,
            8.9411765e-01,
            1.0000000e00,
            0.0000000e00,
            9.7647059e-01,
            1.0000000e00,
            0.0000000e00,
            1.0000000e00,
            9.4509804e-01,
            0.0000000e00,
            1.0000000e00,
            8.7450980e-01,
            0.0000000e00,
            1.0000000e00,
            7.9215686e-01,
            0.0000000e00,
            1.0000000e00,
            7.0588235e-01,
            0.0000000e00,
            1.0000000e00,
            6.0392157e-01,
            0.0000000e00,
            1.0000000e00,
            4.8235294e-01,
            0.0000000e00,
            1.0000000e00,
            3.1372549e-01,
            0.0000000e00,
            1.0000000e00,
            0.0000000e00,
            1.4901961e-01,
            1.0000000e00,
            0.0000000e00,
            3.3333333e-01,
            1.0000000e00,
            0.0000000e00,
            4.4705882e-01,
            1.0000000e00,
            0.0000000e00,
            5.3725490e-01,
            1.0000000e00,
            0.0000000e00,
            6.1176471e-01,
            9.7647059e-01,
            0.0000000e00,
            6.6666667e-01,
            8.9411765e-01,
            0.0000000e00,
            6.6666667e-01,
            7.9607843e-01,
            0.0000000e00,
            6.3921569e-01,
            6.9019608e-01,
            0.0000000e00,
            5.9215686e-01,
            5.6470588e-01,
            0.0000000e00,
            5.0980392e-01,
            3.9607843e-01,
            0.0000000e00,
            3.8039216e-01,
        ]
        colors = np.reshape(colors, (-1, 3))
        name = "flexpart_cmap"
    cmap = ListedColormap(colors, name)
    return cmap
