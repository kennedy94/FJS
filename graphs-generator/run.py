problems= ["DAFJS01.txt", "DAFJS02.txt","DAFJS03.txt",
           "DAFJS04.txt","DAFJS05.txt","DAFJS06.txt","DAFJS07.txt","DAFJS08.txt",
           "DAFJS09.txt","DAFJS10.txt","DAFJS11.txt","DAFJS12.txt","DAFJS13.txt",
           "DAFJS14.txt","DAFJS15.txt","YFJS01.txt", "YFJS02.txt","YFJS03.txt","YFJS04.txt",
           "YFJS05.txt","YFJS06.txt","YFJS07.txt","YFJS08.txt","YFJS09.txt","YFJS10.txt"]


j=1

#m=0.5
#p=0.5
itMax = 1000
delta = 0.99
f=0.0001

import os
import numpy as np
import subprocess
import multiprocessing


i=1
for alpha in np.arange(0.0,1.1,0.2):
    for inst in problems:
        for m in np.arange(0.1,1.0,0.1):
        #for itMax in np.arange(1100,2001,100):
            for p in np.arange(0.1,1.0,0.1):
            #for delta in [0.9, 0.95, 0.99]:
                command = "./FJS"+" "+inst+" SASemBuscaItMax1000Delta099VariandoMeP.csv"+" 600 "+str(-alpha)+" None None Best SA "+str(itMax)+" "+str(m)+" "+str(p)+" "+str(f)+" "+str(delta)+" & \n"
                if i % 10 == 0:
                    command += "wait \n"
                
                with open ('runSASemBuscaItMax1000Delta099VariandoMeP.sh', 'a') as rsh:
                    rsh.write(command)
                i += 1
