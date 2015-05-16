import numpy as np
import matplotlib.pyplot as plt
import pyfits
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy import interpolate, signal
from scipy.interpolate import InterpolatedUnivariateSpline as ius

#### COMBINED K FLUX
combflux = np.loadtxt('combfluxK.txt')
wvdata = pyfits.getdata('SKY_SDCK_20141009_0100.wvlsol_v1.fits')

key=9 # for 2.2 microns
combflux = combflux[key]
wvdata = wvdata[key]
gv = np.isfinite(combflux)
wvdata,combflux=wvdata[gv],combflux[gv]

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print 'x = %f, y = %f'%(ix, iy)
    sx=find_nearest(wvdata, ix)
    ind = np.where(wvdata == sx)
    sy=combflux[ind]
#    l.set_xdata(ix); l.set_ydata(iy)
#    l.set_xdata(sx); l.set_ydata(sy)
    l.set_xdata([ix,sx]); l.set_ydata([iy,sy])
    plt.draw()

    # assign global variable to access outside of function
    global coords
    coords.append((ix, iy))
    return

def fit_spline(coords,wvln,flux):
    sx,sy=np.array([]),np.array([])
    for i in range(len(coords)):
        if coords[i][0] < 1.1: #hacky bug fix; make boolean?
            sx=np.append(sx,np.nan)
            sy=np.append(sy,np.nan)
        else:
            sx=np.append(sx,find_nearest(wvln, coords[i][0]))
            ind = np.where(wvln == sx[i])
            sy=np.append(sy,flux[ind])
    gv = np.isfinite(sy)
    sx,sy=sx[gv],sy[gv]
    print sx, sy
#    splfit=interpolate.splrep(sx,sy,s=0.7) #default deg=3
#    cont = interpolate.splev(wvln,splfit)
    splfit = ius(sx, sy, k=3)
    cont = splfit(wvln)
    return cont, sx, sy

coords = []
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wvdata,combflux,'k-')
l,=ax.plot(wvdata[0],combflux[0],'bo--',markersize=4)
ax.set_xlim([wvdata[0],wvdata[-1]])

# Button to Exit Out and Fit Spline
axcolor = 'lightgoldenrodyellow'
end = plt.axes([0.8, 0.02, 0.15, 0.04])
button = Button(end, 'FIT SPLINE', color=axcolor, hovercolor='0.975')
def end(event):
    fig.canvas.mpl_disconnect(cid)
    plt.close()
    del coords[-1] #hacky bug fix; removes last click
button.on_clicked(end)

# Radio Button to Add/Remove Coordinates for Fitting
#'''
rax = plt.axes([0.025, 0.01, 0.15, 0.07], axisbg=axcolor)
radio = RadioButtons(rax, ('ADD', 'REMOVE'), active=0)
def addrmfunc(action):
    ### commands that enables add/remove points
    # if Add: enable appending coorinates
    # elif Remove: enable removing coorinates (pick nearest coordinate; remove)
    plt.draw()
radio.on_clicked(addrmfunc)
#'''

# Call click function
cid = fig.canvas.mpl_connect('button_press_event', onclick)
ax.set_xlabel(r'$\lambda\ [\mu \rm m]$',size=18)
plt.tight_layout()
plt.show()

cont, xp, yp = fit_spline(coords, wvdata, combflux)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot([wvdata[0],wvdata[-1]],[1,1],color='0.7')
ax.plot(xp,yp,'ro')
ax.plot(wvdata,combflux)
ax.plot(wvdata,cont,'k--',lw=2)
ax.plot(wvdata,combflux/cont)
ax.set_xlim([wvdata[0],wvdata[-1]])
ax.set_xlabel(r'$\lambda\ [\mu \rm m]$',size=18)
plt.tight_layout()
plt.show()
