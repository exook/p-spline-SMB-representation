from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from os import sys, path, remove
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinter
import library
from scipy.stats import chisquare
import pspline8

def fit_atlas(alpha):
    atlas_x,atlas_y=library.atlas()
    (x,y,yd3,atlas_pspline)=library.create_pspline(atlas_x,atlas_y,alpha=alpha)
    chi_atlas=library.chisquare(y,yd3)
    print("chi atlas: ",chi_atlas)
    print('-----------------------------------------------------------')
    return library.makePrettyPlots(x, y, yd3, title="ATLAS, alpha:"
                                   +str(alpha)+", chi:"+str(chi_atlas))


def average_alpha():
    number_of_events=library.number_of_files_in_folder()
    #number_of_events=15
    print("number of runs:", number_of_events)
    list_of_chi=[]
    #alpha=500.625
    alpha=1#1.078125#20:0.875
    mean=0
    factor=1
    operator=''
    while True:
        list_of_chi=[]
        for iteration in range(1,number_of_events+1):
            (x_madgraph,y_madgraph)=library.MadGraph(iteration)
            a=library.create_pspline(x_madgraph,y_madgraph,alpha)
            (x_madgraph,y_madgraph,yd3_madgraph,madgraph_pspline)=a
            list_of_chi.append(library.chisquare(y_madgraph,yd3_madgraph))

        mean=sum(list_of_chi)/number_of_events
        print('Average Chi: ',mean)
        if mean<0.999:
            if operator=="SUBTRACTING":
                factor=factor/2
            operator="ADDING"
            alpha +=0.5*factor
        elif mean>1.001:
            if operator=="ADDING":
                factor=factor/2
            operator="SUBTRACTING"
            alpha =alpha-0.5*factor
        else:
            print('-----------------------------------------------------------')
            print("Best alpha is:",alpha,"\nIt gives a average chi of: ",mean)
            print('-----------------------------------------------------------')
            break
        print("Current alpha: ", alpha)

    atlas_x,atlas_y=library.atlas()
    (x,y,yd3,atlas_pspline)=library.create_pspline(atlas_x,atlas_y,alpha=alpha)
    chi_atlas=library.chisquare(y,yd3)
    
    print("chi atlas: ",chi_atlas)
    print('-----------------------------------------------------------')
    #plta=library.makePrettyPlots(x, y, yd3, title="ATLAS, alpha:"+str(alpha)
    #                             +", chi:"+str(chi_atlas), ymax = 2e5)
    return alpha,list_of_chi

def alpha_spectrum(start,stop,list_of_alpha):
    number_of_events=library.number_of_files_in_folder()
    print("number of runs:", number_of_events)
    alpha=0.5
    factor=2
    operator=''
    for iteration in range(start,stop+1):
        alpha=0.1
        factor=2
        list_of_alpha.append(find_alpha(library.MadGraph(iteration),
                                        alpha,factor,operator))
    return list_of_alpha

def find_alpha(tuuple,alpha,factor,operator):
    while True:
        x,y=tuuple
        (x,y,yd3,pspline)=library.create_pspline(x,y,alpha)
        this_chi=(library.chisquare(y,yd3))
        #print('Run: ',iteration,' Chi: ',this_chi,' alpha: ',alpha)
        if 0.999<this_chi<1.001:
            break
        elif this_chi<0.999:
            if operator=="SUBTRACTING":
                factor=factor/2
            operator="ADDING"
            alpha +=0.5*factor
        elif this_chi>1.001:
            if operator=="ADDING":
                factor=factor/2
            operator="SUBTRACTING"
            alpha =alpha-0.5*factor
            if alpha<0:
                while alpha<0:
                    alpha=alpha+1.5
                    factor=factor/2
    return alpha    

def plot_alpha_spectrum(list_of_alpha,atlas_alpha):
    bin_list = np.linspace(0, 10, 51)
    plt.rc("figure", facecolor="white")
    plt.hist(list_of_alpha, bins=bin_list)
    # plt.hist passes it's arguments to np.histogram
    plt.title("Alpha distribution\nOptimal ATLAS " r'$\alpha$''='
              +str(round(atlas_alpha,2)))
    plt.axis([0, 10, 0, 15])
    plt.ylabel("Counts")
    plt.xlabel(r"$\alpha$")
    extraticks=[atlas_alpha]
    plt.xticks(list(plt.xticks()[0]) + extraticks)
    print(sum(list_of_alpha)/len(list_of_alpha))
    return plt

def plot_chi_spectrum(list_of_alpha,atlas_alpha):
    plt.rc("figure", facecolor="white")
    bin_list = np.linspace(0, 2, 21)
    plt.hist(list_of_alpha, bins=bin_list)
    # plt.hist passes it's arguments to np.histogram
    plt.title(r"$\chi^2$"+" Distribution")
    plt.ylabel("Counts")
    plt.xlabel(r"$\chi^2$")
    
    plt.axis([0, 2, 0, 30])
    extraticks=[atlas_alpha]
    plt.xticks(list(plt.xticks()[0]) + extraticks)
    return plt


def makePrettyPlots(xs, ys, bkgs, title,ymin=1e-1, ymax = 1e5,xmin=1100):
    #This function is credited to Meghan Frate
    f, (ax1, ax2) = plt.subplots(2, sharex=True,
                                 figsize=(12,12),
                                 gridspec_kw = {'height_ratios':[3, 1]})
    f.suptitle(title, fontsize=30)
    dataPlot = ax1.errorbar(xs, ys, marker='o', ls='None',
                            yerr = np.sqrt(ys), c='black', markersize=7, label="Data")
    bkg1Plot, = ax1.plot(xs, bkgs, color='g', linewidth=3.0, label="P-Spline")

    ax1.legend()
    ax1.set_ylabel('Number of Events', fontsize=20)
    ax1.set_xlabel('m (GeV)', fontsize=20)
    #ax1.set_yscale('log', nonposy="clip")
    #ax1.set_xscale('log')
    ax1.set_xlim([xmin, 7500])
    ax1.set_ylim([ymin, ymax])
    ax1.tick_params(axis='y', labelsize=20)

    diff=(ys-bkgs)/bkgs
    ax2.plot(xs, diff, marker='o', ls='None', c='black', markersize=7,)
    ax2.axhline(0, color='black', lw=1)
    ax2.set_ylabel(r'$\frac{y_{Data}-y_{spline}}{y_{spline}}$', fontsize=20)
    ax2.set_xlabel('Invariant Mass (GeV)', fontsize=20)
    #ax2.set_xscale('log')
    ax2.tick_params(axis='y', labelsize=20)
    ax2.set_xlim([xmin, 7500])
    ax2.set_ylim([-1.4, 1.4])


    ax2.set_xticks([xmin, 2000 ,3000, 4000, 5000, 6000, 7000])
    minor_ticks = np.arange(xmin, 7500, 100)
    ax2.set_xticks(minor_ticks, minor=True)  
    labels = [str(xmin),str(xmin+900),str(xmin+1900),
              str(xmin+2900),str(xmin+3900),str(xmin+4900),str(xmin+5900)]
    ax2.set_xticklabels(labels)
    
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)   
  
    #plt.show()
    return plt

def atlas_no():
    read_atlas_x=open("xval_ATLAS.txt","r")
    x=[]
    for line in read_atlas_x:
        x.append(float(line[:-1]))
    read_atlas_x.close()
    x=np.array(x)

    read_atlas_y=open("yval_ATLAS.txt","r")
    y=[]
    for line in read_atlas_y:
        y.append(float(line[:-1]))
    read_atlas_y.close()
    y=np.array(y)

    return x,y

def atlas():
    ATLAS_hist=open('ATLAS_hist.txt','r')
    data=[]
    for line in ATLAS_hist:
            data.append(float(line[:-1]))
            
    read_bins=open("original_bins.txt","r")
    bins=[]
    for line in read_bins:
        bins.append(float(line[:-1]))
    
    hist=np.histogram(data,bins=bins, range=None,
                      normed=False, weights=None, density=None)
    y=hist[0]
    read_bins.close()
    ATLAS_hist.close()
    return (hist[1][:-1],y)

def create_pspline(x,y,alpha):
# This function is heavily based on a file that is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# For this function:
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
    a=splinter.BSplineBuilder(x, y, smoothing=splinter.BSplineBuilder.Smoothing.PSPLINE,
                              alpha=alpha).build()
    pspline=a
    n = len(x)
    xd = [0.]*n
    yd3 = [0.]*n
    for i,xd in enumerate(x):
        yd3[i] = pspline.eval(xd)[0]

    yd3=np.array(yd3)
    return (x,y,yd3,pspline)

def MadGraph(iteration):#with bin width division

    p=Path("/media/alexander/INTENSO/gp_events/run_"
           +str(iteration)+"/invariant_masses_final.txt")
    invariant_masses=p.open("r")
    
    data=[]
    for line in invariant_masses:
            data.append(float(line[:-1]))
    read_bins=open("original_bins.txt","r")
    
    bins=[]
    for line in read_bins:
        bins.append(float(line[:-1]))
    
    hist=np.histogram(data,bins=bins, range=None,
                      normed=False, weights=None, density=None)

    y=hist[0]
    read_bins.close()
    invariant_masses.close()
    return (hist[1][:-1],y)

def chisquare(y,yd3):
    chisquare=(sum(((y-yd3)**2)/yd3))/(len(y)-1)
    return chisquare

def number_of_files_in_folder():
    p=Path("/media/alexander/INTENSO/gp_events")
    return len(list(p.iterdir()))

def plot_all_mg(number_of_events,alpha):
    number_of_events=library.number_of_files_in_folder()
    for iteration in range(1,number_of_events+1):
            (x_madgraph,y_madgraph)=library.MadGraph(iteration)
            a=library.create_pspline(x_madgraph,y_madgraph,alpha)
            (x_madgraph,y_madgraph,yd3_madgraph,madgraph_pspline)=a
            pltb=library.makePrettyPlots(x_madgraph, y_madgraph, yd3_madgraph,
                            title="MadGraph run"+str(iteration)
                            +", alpha:"+str(alpha)+", chi:"
                            +str(library.chisquare(y_madgraph,yd3_madgraph)),
                            ymax = 2e5) 
            pltb.show()

def plot_one_mg(run_number,alpha):
    (x_madgraph,y_madgraph)=library.MadGraph(run_number)
    a=library.create_pspline(x_madgraph,y_madgraph,alpha)
    (x_madgraph,y_madgraph,yd3_madgraph,madgraph_pspline)=a
    pltb=library.makePrettyPlots(x_madgraph, y_madgraph, yd3_madgraph,
                    title="MadGraph run"+str(run_number)+", alpha:"+str(alpha)
                    +", chi:"+str(library.chisquare(y_madgraph,yd3_madgraph)),
                    ymax = 2e5)
    return pltb

    
