import os
import numpy as np
import itertools as itert
from daskms.experimental.zarr import xds_from_zarr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams
import matplotlib.dates as mdates
import pytz
from pyrap.tables import table
from astropy.time import Time
from scipy.constants import c as c_speed
from scipy.interpolate import interp1d
from pytz import timezone
from datetime import datetime
plt.rcParams["figure.figsize"] = [12, 9]
plt.rcParams["figure.dpi"] = 200



SMALL_SIZE = 20
MEDIUM_SIZE = 16
BIGGER_SIZE = 24


plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#**************************************************
def plot_freq_vs_time(path_out, arr0, f_dict, scann, ref_ant, label0="K", corr=0, ftick0=100, ftick1=200, \
    zoomed=False, arr2=None):
    """
    Imshow amp/ phase of gains.

    """

    #Get the antenna info.
    antennas = arr0[scann].antenna


    #Get the xaxis ticks
    time = arr0[scann]["gain_time"]
    ntime = time.shape[0]


    #Time converted to our calendar (UTC)
    utc_times = arr_obs_time(time)


    #Change time to SAST.
    # sast_zone = pytz.timezone("Africa/Johannesburg")
    gmt2_zone = pytz.timezone("Europe/Amsterdam")

    if ntime != 1:
        #Reduced number of ticks for x axis.
        tick_num0 = time.shape[0]// 2 #thus, ideally only showing three ticks on the x axis.
        tick_num0 += 1
        reduced_ticks0 = np.arange(0, len(utc_times), tick_num0) #show every tick_num0 ticks
        labels0 = np.array([str(utc_times[i]) for i in reduced_ticks0])

        #Only show hours and minutes
        reduced_ticklabels0 = [ \
            (datetime.strptime(label, "%Y-%m-%d %H:%M:%S.%f").astimezone(gmt2_zone)).strftime("%H:%M") for label in labels0
            ]
        
        #x axis title
        # xtitle = "Time (SAST)"
        xtitle = "Time (GMT+2)"
    else:
        xtitle = "Time"


    #Get the y axis.
    freq0 = arr0[scann]["gain_freq"].values
    freq = 1e-6 * freq0 #in MHz
    nchan = len(freq)

    if nchan != 1:
        f11 = ftick0
        f22 = ftick1
        f1 = (f11 - freq[0])/ (freq[-1] - freq[0])
        f2 = (f22 - freq[0])/ (freq[-1] - freq[0])

        reduced_ticks1 = np.array([f1*nchan, f2*nchan])
        reduced_ticklabels1 = np.array([f11, f22])

        #y axis title
        ytitle = "Frequency (MHz)"
    else:
        ytitle = "Frequency"



    #looking at all antennas
    n_ant = len(antennas)
    nrows = int(np.sqrt(n_ant))
    ncols = n_ant//nrows
    if nrows*ncols != n_ant:
        ncols += 1
    

    #Choose the correlation
    if corr == 0:
        label1 = "XX"
    elif corr == 3:
        label1 = "YY"

    #Compute referenced gains.
    gains = make_referenced_gains_qc(arr0, scann, ref_ant=ref_ant)
    

    #Working on the gain flags
    if arr2 is not None:
        #define the difference of the phases
        gains2 = make_referenced_gains_qc(arr2, scann, ref_ant=ref_ant)
        
        arr1 = f_dict["func"](gains.flatten()*np.conjugate(gains2.flatten().T))
        arr1 = arr1.reshape(gains.shape)
        #Instead of using or, use logical_or for sorting through different arrays
        sel = np.where(np.logical_or(arr0[scann].gain_flags.values == 1, arr2[scann].gain_flags.values == 1))

    else:
        arr1 = f_dict["func"](gains)
        ##Get the gains_flags.
        sel = np.where(arr0[scann].gain_flags.values == 1)

    
    #Zoomed in colourbar.
    if zoomed:

        vmax = 0.5
        vmin = -0.5
        label2 = "_zoomed"
    else:
        vmax = 1.01 * arr1.max()
        vmin = 0.99 * arr1.min()
        label2 = ""


    ind = 0
    fig, axs = plt.subplots(nrows, ncols, figsize=(24, 18), sharex=True, sharey=True, constrained_layout=True)
    for i1 in range(nrows):
        for i2 in range(ncols):
            if ind <= n_ant-1:
                if f_dict["func"] == np.angle:
                    cmap = matplotlib.colormaps["hsv"]
                else:
                    cmap = "viridis"

                arr1[sel] = np.nan
                im = axs[i1, i2].imshow(arr1[:, :, ind, 0, corr].T, cmap=cmap, interpolation="None", \
                    aspect="auto", origin="lower", vmin=vmin, vmax=vmax)

                if ntime != 1:
                    axs[i1, i2].set_xticks(reduced_ticks0)
                    axs[i1, i2].set_xticklabels(reduced_ticklabels0)


                if nchan != 1:
                    axs[i1, i2].set_yticks(reduced_ticks1)
                    axs[i1, i2].set_yticklabels(reduced_ticklabels1)


                #Label each subplot with the antenna indices.
                axs[i1, i2].text(0.02, 0.92, "%s"%antennas.values[ind], color="k", transform=axs[i1, i2].transAxes)

            ind += 1

    fig.supxlabel(xtitle)
    fig.supylabel(ytitle)
    fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.6)
    fig.savefig(path_out+"{}_freq_time_{}{}{}_scan{}.pdf".format(label0, f_dict["label1"], label1, label2, scann), \
        bbox_inches="tight")
    plt.close()



def make_referenced_gains_qc(arr0, scann, ref_ant):
    """
    Return a numpy array containing referenced gains with 
    respect to the reference antenna.

    """

    arr0 = arr0[scann]
    gains0 = arr0.gains.values
    antennas = arr0.antenna
    
    #The number of antennas.
    n_ant = len(antennas)
    
    for p in range(n_ant):
        gains0[:, :, p] *= np.exp(-1.0j*np.angle(gains0[:, :, ref_ant]))    

    return gains0


def arr_obs_time(time):
    """
    Convert Modified Julian Dates (seconds) to UTC time zone.

    """
    #The values from QC are in MJD in seconds.
    #Hence, first of all, convert the values into days.
    time_mjd = time/(60*60*24)
    time_jd = time_mjd + 2400000.5 #change to Julian date

    #Change to our calendar dates.
    time_calendar = Time(time_jd, format="jd")
    
    #Time in utc time zone
    return time_calendar.utc.iso