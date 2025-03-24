# PlottiQal
PlottiQal is not an official [QuartiCal](https://github.com/ratt-ru/QuartiCal) project. It contains custom scripts for visualising and analysing QuartiCal outputs.

Usage: main.py [options]

Options:
  -h, --help            show this help message and exit
  --plotop=PLOTOP       Plotting options: gains, params or all
  --sol1=SOL1           QuartiCal solutions (xarray Dataset)
  --sol2=SOL2           Difference with solutions sol
  --ref_ant=REF_ANT     Reference antenna
  --corr=CORR           Specify correlation axis
  --scann=SCANN         Specify scan number
  --label0=LABEL0       Figure label
  --plotpar=PLOTPAR     Specify parameter to be plotted (delay or TEC)
  --par_ind=PAR_IND     Specify parameter index
  --ref_nparr=REF_NPARR
                        Use to compare parameter plots (same time axis)
  --ftick0=FTICK0       Specify frequency axis tick 1
  --ftick1=FTICK1       Specify frequency axis tick 2
  --zoomed=ZOOMED       Additional zoomed in plot on specific colourbar
  -o OPDIR, --opdir=OPDIR
                        Output folder to store plots
