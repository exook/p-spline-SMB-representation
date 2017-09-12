from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from os import sys, path, remove
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinter
from scipy.stats import chisquare
import library

def test():
    alpha=0.8187500000000001
    atlas_x,atlas_y=library.atlas()
    (x,y,yd3,atlas_pspline)=library.create_pspline(atlas_x,atlas_y,alpha=alpha)
    chi_atlas=library.chisquare(y,yd3)
    print("chi atlas: ",chi_atlas)
    print('-----------------------------------------------------------')
    plta=library.makePrettyPlots(x, y, yd3, title="ATLAS data fitted with "
    +r"$\alpha=$"+str(round(alpha,2))+"\n"+r'$\chi^2=$'+str(round(chi_atlas,2)))
    list_of_chi_final=[]
    iteration=18
    (x_madgraph,y_madgraph)=library.MadGraph(iteration)
    a=library.create_pspline(x_madgraph,y_madgraph,alpha)
    (x_madgraph,y_madgraph,yd3_madgraph,madgraph_pspline)=a
    chi_madgraph=library.chisquare(y_madgraph,yd3_madgraph)
    plta.show()

def run():
    alpha,list_of_chi=library.average_alpha()
    pltc=library.plot_chi_spectrum(list_of_chi,1)
    pltc.show()
    #pltc=library.plot_alpha_spectrum(list_of_chi,1)
    #pltc.show()
    list_of_alpha=[]
    list_of_alpha=library.alpha_spectrum(1,100,list_of_alpha)#inclusive
    atlas_alpha=library.find_alpha(library.atlas(),0.1,2,'')
    print('Atlas has a optimal alpha of: ',atlas_alpha)    
    pltb=library.plot_alpha_spectrum(list_of_alpha,atlas_alpha)
    pltb.show()
    pltd=library.fit_atlas(atlas_alpha)
    pltd.show()

if __name__=="__main__":
    run()
    #test()


    
    



