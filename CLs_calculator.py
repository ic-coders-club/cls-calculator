#! /usr/bin/env python
import numpy as np
from scipy.stats import poisson
from scipy.stats import norm

#rename function
Poiss=poisson.pmf

#set some parameters
ntoys=10000

#get data from text file
a=np.loadtxt('data.txt',comments='#',ndmin=2)
bg=a[:,0]
bg_unc=a[:,1]
sig=a[:,2]
sig_unc=a[:,3]
data=a[:,4]
nbins=len(a)

#the test statistic
def t(x,bg,sig):
    q=Poiss(x,bg+sig)/Poiss(x,bg)
    return -2*np.sum(np.log(q),axis=1)

#generate the ntoys*nbins means of the poisson distribution 
bg_for_toys=norm.rvs(loc=bg,scale=bg*bg_unc,size=(ntoys,nbins))
sig_for_toys=norm.rvs(loc=sig,scale=sig*sig_unc,size=(ntoys,nbins))
#throw one poisson per mean 
bg_toys=poisson.rvs(bg_for_toys) 
sigbg_toys=poisson.rvs(bg_for_toys+sig_for_toys) 
#calculate the test statistic
t_bg_toys=t(bg_toys,bg,sig)
t_sigbg_toys=t(sigbg_toys,bg,sig)
t_data=t(np.array(data,ndmin=2),bg,sig)
#calculate CLS
nbg=float(len(np.where(t_bg_toys>t_data)[0]))
nsigbg=float(len(np.where(t_sigbg_toys>t_data)[0]))
print('t={}'.format(t_data[0]))
cls=nsigbg/nbg
print('cls: {}'.format(cls))

