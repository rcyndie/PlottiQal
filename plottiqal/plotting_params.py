import os
import numpy as np
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
# import juliandate as jd
from datetime import datetime
plt.rcParams["figure.figsize"] = [12, 9]
plt.rcParams["figure.dpi"] = 200

# #other thesis plots
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

#**********************************************************************************************
def plot_params(path_out, arr0, scann, ref_ant, label0="K", plotpar="delay", par_ind=1, ref_nparr=None):
    """
    Plot params from gains.params.
    (par_ind = param index)

    """

    params = arr0[scann].params.values
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
        reduced_ticks0 = np.arange(0, ntime, tick_num0) #show every tick_num0 ticks
        labels0 = np.array([str(utc_times[i]) for i in reduced_ticks0])

        #Only show hours and minutes
        reduced_ticklabels0 = [ \
            (datetime.strptime(label, "%Y-%m-%d %H:%M:%S.%f").astimezone(gmt2_zone)).strftime("%H:%M") for label in labels0
            ]
    xtitle = "Time (min)"


    ##Get the gains_flags.
    sel = np.where(arr0[scann].param_flags.values == 1)


    #looking at all antennas
    n_ant = len(antennas)
    nrows = int(np.sqrt(n_ant))
    ncols = n_ant//nrows
    if nrows*ncols != n_ant:
        ncols += 1


    if plotpar == "TEC":
        rescaling_const = (-c_speed/40.308)*1e-16
        ytitle1 = "dTEC (TECU)"
        colour_line = "blue"
        legend0 = r"$\text{T}_\text{QC}$"
    elif plotpar == "delay":
        rescaling_const = 1.
        ytitle1 = "Delay (s)"
        colour_line = "k"
        legend0 = r"$\text{K}_\text{QC}$"


    ##ind is indexing across antennas
    ind = 0
    fig, axs = plt.subplots(nrows, ncols, figsize=(24, 18), sharex=True, sharey=True, constrained_layout=True)
    for i1 in range(nrows):
        for i2 in range(ncols):
            if ind <= n_ant-1:
                params[sel] = np.nan
                param = rescaling_const * (params[:, 0, ind, 0, par_ind]-params[:, 0, ref_ant, 0, par_ind])
                axs[i1, i2].plot(utc_times, param, label=legend0, linewidth=4., alpha=0.5, color=colour_line)

                if ref_nparr is not None:
                    axs[i1, i2].plot(utc_times, ref_nparr[:, ind]-ref_nparr[:, ref_ant], linewidth=1.5, linestyle="--", \
                        color=colour_line, label="DP3")

                axs[i1, i2].axhline(y=0., linestyle=":")


                axs[i1, i2].text(0.02, 0.92, "%s"%antennas.values[ind], color="k", transform=axs[i1, i2].transAxes)

                if ntime != 1:
                    axs[i1, i2].set_xticks(reduced_ticks0)
                    axs[i1, i2].set_xticklabels(reduced_ticklabels0)

                if plotpar == "TEC":
                    # axs[i1, i2].set_ylim(-0.01, 0.01)
                    # axs[i1, i2].set_ylim(-0.5, 0.5)
                    axs[i1, i2].set_ylim(-0.75, 0.75)

                elif plotpar == "delay":
                    # axs[i1, i2].set_ylim(-1e-7, 1e-7)
                    axs[i1, i2].set_ylim(-2.5e-7, 1e-7)

                if ind == ref_ant:
                    axs[i1, i2].legend(loc="best")


            ind += 1

    fig.supxlabel(xtitle)
    fig.supylabel("{}".format(ytitle1))
    fig.savefig(path_out+"{}_params_ind{}_scan{}.pdf".format(label0, par_ind, scann), bbox_inches="tight")
    plt.close()


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