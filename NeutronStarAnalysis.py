
import sys

import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
import astropy
from astropy.io import fits
import regions
import skyfield as sky
#from xpbin import xpbin as _xpbin
import sys
import os
sys.path.insert(0,os.path.abspath(f'/home/jacob/Desktop/HerX-1/02004101'))
import event_l1, event_l2, hk, auxil
import math
#from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
#if PYXSPEC_INSTALLED:
#    import ixpeobssim.evt.xspec_.py as xspec
import pandas as pd
import ixpeobssim as ixpe 
from scipy import optimize as opt

import pickle
from scipy import pi as pi
from scipy.special import jn  as jn
import scipy.interpolate, scipy.misc
import IPython.display as ipd
from scipy import signal
from scipy.interpolate import InterpolatedUnivariateSpline
import xlsxwriter
import openpyxl
from glob import glob

def ImportOrbitalCorrectionData(OBSID,SOURCENAME,SOURCERAD):
    orbdatafile = f'/home/jacob/Desktop/{SOURCENAME}/herx1.fits'
    orbdataopen = astropy.io.fits.open(orbdatafile)
    orbdata = orbdataopen[1].data
    spindata = orbdataopen[2].data
    orbdataopen.close()

    RA = orbdata['RA']
    DEC = orbdata['DEC']
    BINARY = orbdata['BINARY']
    PBINARY = orbdata['PBINARY']
    PBDOT = orbdata['PBDOT']
    BINARYEPOCH = orbdata['BINARYEPOCH']
    AXSINI = orbdata['AXSINI']
    ECCENTRICITY = orbdata['ECCENTRICITY']
    EGRESS = orbdata['EGRESS']
    INGRESS = orbdata['INGRESS']
    OMEGA = orbdata['PERIAPSE']


#def ImportModulationFactor(OBSID,SOURCENAME,SOURCERAD):
    modulationfiles = sorted(glob(f'/home/jacob/Desktop/{SOURCENAME}/RFiles/ModulationFiles/ixpe_d*mfact_02.fits'))
    for i,f in enumerate(modulationfiles):
        detcode = int(f[f.find('ixpe_d')+6])
        print(i,f,detcode)
        with fits.open(modulationfiles[detcode-1]) as hdu:
            if i == 0:
                moddata = hdu[1].data
                moddetector = np.full(len(hdu[1].data),detcode)
                hdu.close()
            else:
                moddata = np.append(moddata,hdu[1].data)
                moddetector = np.append(moddetector,np.full(len(hdu[1].data),detcode))
                hdu.close()
    moddata = np.array([moddata['ENERG_LO'],moddata['ENERG_HI'],moddata['SPECRESP'],moddetector]).T
    moddata = pd.DataFrame(moddata,columns=['ENERG_LO', 'ENERG_HI','SPECRESP','detector'])


#def ImportDataFiles(OBSID,SOURCENAME,SOURCERAD):
    BinaryCorrection = True
    datafiles = sorted(glob(f'/home/jacob/Desktop/HerX-1/{OBSID}/event_l2/*barycorr*.fits*'))
    for i,f in enumerate(datafiles):
        detcode = int(f[f.find('det')+3])
        print(i,f,detcode)
        with fits.open(datafiles[detcode-1]) as hdu:
            if (i ==0):
                MJDREFI=hdu[0].header['MJDREFI']
                MJDREFF=hdu[0].header['MJDREFF']
                LIVETIME = hdu[0].header['LIVETIME']
                events = hdu[1].data
                detector = np.full(len(hdu[1].data),detcode)
                hdu.close()
            else:
                events = np.append(events,hdu[1].data)
                detector = np.concatenate((detector,np.full(len(hdu[1].data),i+1)))
                hdu.close()
#def SetupBackGroundandEnergyBounds(OBSID,SOURCENAME,SOURCERAD):
    coords = pd.read_csv(f'/home/jacob/Desktop/{SOURCENAME}/CoordinatesofNS.txt',delim_whitespace=True)
    CenterX, CenterY= coords['X'].loc[coords['OBS_ID'] == float(OBSID)].mean(),coords['Y'].loc[coords['OBS_ID'] == float(OBSID)].mean()
    print(CenterX,CenterY)
    radius = np.sqrt(((events['X']-CenterX)**2)+((events['Y']-CenterY)**2))

    Bounds = (events['PI']>49) & (events['PI']<=200) & (radius<SOURCERAD)
    Background = (events['PI']>49) & (events['PI']<=200) & (radius<4*SOURCERAD) & (radius>2*SOURCERAD)


    evtback = events[Background]
    events = events[Bounds]
    print(len(events))
    detback = detector[Background]
    detector = detector[Bounds]

#def BinaryCorrection(OBSID,SOURCENAME,SOURCERAD):
    day=3600*24



    TPIOVER2O=orbdata['BINARYEPOCH']-2400000.5
    PORB=orbdata['PBINARY']*day
    PORBDOT=orbdata['PBDOT']
    ECCENTRICITY=orbdata['ECCENTRICITY']
    OMEGA=orbdata['PERIAPSE']
    ASINIOVERC=orbdata['AXSINI']
    NUORB=1/PORB
    NUORBDOT=-1/PORB/PORB*PORBDOT
    met0=events['TIME'].min()##evtlist['TIME'].min()
    MJDSTART=MJDREFF+MJDREFI+met0/day
    print(MJDSTART)
    NU=np.interp(MJDSTART,spindata['BARYTIME'],spindata['FREQUENCY'])
    NUDOT=np.interp(MJDSTART,spindata['BARYTIME'],spindata['FDOT'])
    TREFNU=MJDSTART
    torb=((MJDREFF+MJDREFI)+events['TIME']/day-TPIOVER2O)*day##((MJDREFF+MJDREFI)+evtlist['TIME']/day-TPIOVER2O)*day
    phaseorb=(0.25+torb*(NUORB+torb*NUORBDOT/2))                    
    orbtime=ASINIOVERC*((np.sin(2*np.pi*phaseorb)-np.sin(np.pi/2))-
                    0.5*ECCENTRICITY*(np.sin(4*np.pi*phaseorb-np.radians(OMEGA))-
                                        np.sin(np.pi/2-np.radians(OMEGA))))
    deltat=events['TIME']-orbtime-met0#evtlist['TIME']-orbtime-met0
    freqtst=NU
    if (NUDOT==0) & (MJDSTART < spindata['BARYTIME'].max()):
        _ii=spindata['BARYTIME']>59000
        nufunk=InterpolatedUnivariateSpline(spindata['BARYTIME'][_ii],spindata['FREQUENCY'][_ii],k=2)
        nudd=nufunk.derivatives(MJDSTART)
        NU=nudd[0]
        NUDOT=nudd[1]/day
        NUDDOT=nudd[2]/day/day
    elif (NUDOT==0) & (MJDSTART > spindata['BARYTIME'].max()):
        NU = spindata['FREQUENCY'][-1]
        #NUDOT = spindata['FDOT'][-1]
        NUDDOT = 0
    else:
        NUDDOT=0

    #Adjust NU + (x) until the 2dhist looks aligned

    adjustment = 0#1e-6#0#5e-5#1e-6#1e-7
    phase=np.mod(deltat*((NU+adjustment)+0.5*deltat*(NUDOT*deltat*NUDDOT*0.3333333)),1)


    df = np.array([events['TIME'],events['PI'],events['U'],events['Q'],events['W_MOM'],events['X'],events['Y']]).T
    #df = df.byteswap().newbyteorder()

    df = pd.DataFrame(df,columns=['Barycentered Times','CH','U','Q','W_mom','X','Y'])
    df['Time'] = deltat
    df['Phase'] = phase
    df['detector'] = detector
    df['OrbitalCorrectedTime'] = df['Time']

    bindata=plt.hist2d(df['Phase'],df['Time'],bins=100)
    plt.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/BeginningPeriodCorrection.pdf')
    #for i in range(1,4):
    #    print(i)
    #    df.loc[df['detector']==i].to_csv(f'/home/jacob/Desktop/{OBSID}dataframeDetector{i}.csv',sep=',',header=True,mode='w',index=False,columns=['Barycentered Times','OrbitalCorrectedTime','Phase','X','Y'])


#def PeriodCorrection(OBSID,SOURCENAME,SOURCERAD):
    def FindPeriod(adjustment):
        phase=np.mod(deltat*((NU+adjustment)+0.5*deltat*(NUDOT*deltat*NUDDOT*0.3333333)),1)
        PulseP = np.histogram(phase,bins=100)
        PF = (PulseP[0].max()-PulseP[0].min())/(PulseP[0].max()+PulseP[0].min())
        return 1/PF
    bnds = [(-1e-4,1e-4)]
    initialguess=[0]
    res = opt.minimize(FindPeriod,initialguess,method='nelder-mead',bounds = bnds)
    NU += res.x
    print(NU)
    phase=np.mod(deltat*((NU)+0.5*deltat*(NUDOT*deltat*NUDDOT*0.3333333)),1)
    df['Phase'] = phase
    bindata=plt.hist2d(df['Phase'],df['Time'],bins=100)
    plt.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/FinalPeriodCorrection.pdf')


#def PulseProfile(OBSID,SOURCENAME,SOURCERAD):
    PulseP = np.histogram(phase,bins=100)
    print((PulseP[0].max()-PulseP[0].min())/(PulseP[0].max()+PulseP[0].min()))
    Counts = PulseP[0]
    P = []

    CenterofBins = (PulseP[1][:-1] + PulseP[1][1:])/2
    CenterofBins_ = np.append(CenterofBins,CenterofBins+1)
    PulseP_ = np.append(PulseP[0],PulseP[0])

    plt.figure()
    plt.plot(CenterofBins_,PulseP_/LIVETIME)
    plt.xlabel('Phase')
    plt.ylabel('COUNTS/S')
    plt.title(f'{OBSID} Pulse Profile')
    #plt.ylim(0,0.06)
    fig1 = plt.gcf()
    plt.show
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/PulseProfile.pdf')

#def LightCurve(OBSID,SOURCENAME,SOURCERAD):
    binwidth = 500
    lightcurve = plt.hist(df['Time'], bins = range(0,int(df['Time'].max()),binwidth) )
    Counts = lightcurve[0]
    CenterofLCurveBins = (lightcurve[1][:-1] + lightcurve[1][1:])/2

    plt.figure()
    plt.scatter(CenterofLCurveBins,Counts/binwidth,s=3)
    plt.xlabel('Time [s]')
    plt.ylabel('COUNTS/S')
    plt.title(f'{OBSID} Light Curve')
    plt.ylim(0.01,10)
    plt.yscale('log')
    fig1 = plt.gcf()
    plt.show
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/LightCurve.pdf')

    #NOTE for later, the lowered bits are caused by off times within the 250s bins, use the ontime columns to rectify

#def ImportARFMRF(OBSID,SOURCENAME,SOURCERAD):
    ARFfiles = sorted(glob(f'/home/jacob/Desktop/HerX-1/RFiles/*.arf'))

    for i,f in enumerate(ARFfiles):
        detcode = int(f[f.find('ixpe_d')+6])
        print(i,f,detcode)
        with fits.open(ARFfiles[detcode-1]) as hdu:
            if (i ==0):
                Response = hdu[1].data
                detector = np.full(len(hdu[1].data),detcode)
            else:
                Response = np.append(Response,hdu[1].data)
                detector = np.concatenate((detector,np.full(len(hdu[1].data),i+1)))

    Response = np.array([Response['ENERG_LO'].flatten(),Response['ENERG_HI'].flatten(),Response['SPECRESP'].flatten()]).T
    #Response = Response.byteswap().newbyteorder()
    ARF = pd.DataFrame(Response,columns=['E_low','E_hi','Resp'])
    ARF['detector'] = detector

    RMFfiles = sorted(glob(f'/home/jacob/Desktop/{SOURCENAME}/RFiles/*.rmf'))

    for i,f in enumerate(RMFfiles):
        detcode = int(f[f.find('ixpe_d')+6])
        print(i,f,detcode)
        with fits.open(RMFfiles[detcode-1]) as hdu:
            if (i ==0):
                Response = hdu[2].data
                detector = np.full(len(hdu[2].data),detcode)
            else:
                Response = np.append(Response,hdu[2].data)
                detector = np.concatenate((detector,np.full(len(hdu[2].data),i+1)))
    Response = np.array([Response['E_MIN'].flatten(),Response['E_MAX'].flatten(),Response['CHANNEL'].flatten()]).T
    #Response = Response.byteswap().newbyteorder()
    RMF = pd.DataFrame(Response,columns=['E_low','E_hi','CH'])
    RMF['E_AVG'] = (RMF['E_low'] + RMF['E_hi'])/2
    RMF['detector'] = detector

#def ChanneltoEnergyandWeights(OBSID,SOURCENAME,SOURCERAD):
    RMF = RMF.merge(moddata, left_on=['E_low','E_hi','detector'], right_on=['ENERG_LO','ENERG_HI','detector'], how='left' )
    RMF = RMF.drop(columns=['ENERG_LO','ENERG_HI'])


    df = df.merge(RMF,left_on=['CH','detector'],right_on=['CH','detector'],how='left')
    df = df.drop(columns=['E_low','E_hi'])

    df['U*W_mom'] = df['U']* df['W_mom']
    df['Q*W_mom'] = df['Q']* df['W_mom']

#def ISpectra(OBSID,SOURCENAME,SOURCERAD):
    binwidth = 0.04 #The channel energy width is 0.04
    ISpectra = np.histogram(df['E_AVG'], bins = np.arange(2,8,binwidth) )
    Counts = ISpectra[0]
    CenterofISpectraBins = (ISpectra[1][:-1] + ISpectra[1][1:])/2

    plt.figure()
    plt.scatter(CenterofISpectraBins,Counts/LIVETIME/binwidth,s=3)
    plt.xlabel('Energy [KeV]')
    plt.ylabel('COUNTS/S/KeV')
    plt.title(f'{OBSID} I Spectra')
    plt.ylim(0.01,10)
    plt.yscale('log')
    fig1 = plt.gcf()
    plt.show()
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/CombinedISpectra.pdf')

    binwidth = 0.04 #The channel energy width is 0.04
    #ISpectra = plt.hist(df['E_AVG'], bins = np.arange(2,8,binwidth) )
    plt.figure()
    for det in range(1,3+1): #plus 1 so that we get all 3 detectors

        ISpectra = np.histogram(df['E_AVG'].loc[df['detector'] == det], bins = np.arange(2,8,binwidth) )
        Counts = ISpectra[0]
        CenterofISpectraBins = (ISpectra[1][:-1] + ISpectra[1][1:])/2


        plt.scatter(CenterofISpectraBins,Counts/LIVETIME/binwidth,s=3,label=f'Detector: {det}')
    plt.xlabel('Energy [KeV]')
    plt.ylabel('COUNTS/S/KeV')
    plt.title(f'{OBSID} I Spectra')
    #plt.ylim(0.01,10)
    plt.yscale('log')
    plt.legend()
    fig1 = plt.gcf()
    plt.show
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/ISpectra.pdf')

    binwidth = 0.04 #The channel energy width is 0.04
    #ISpectra = plt.hist(df['E_AVG'], bins = np.arange(2,8,binwidth) )
    plt.figure()
    for det in range(1,3+1): #plus 1 so that we get all 3 detectors

        ISpectra = np.histogram(df['E_AVG'].loc[df['detector'] == det], bins = np.arange(2,8,binwidth) )
        Counts = []
        for i in np.arange(2,8-binwidth,binwidth):
            C = df['W_mom'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].sum()
            Counts = np.append(Counts,C)
        CenterofISpectraBins = (ISpectra[1][:-1] + ISpectra[1][1:])/2


        plt.scatter(CenterofISpectraBins,Counts/LIVETIME/binwidth,s=3,label=f'Detector: {det}')
    plt.xlabel('Energy [KeV]')
    plt.ylabel('COUNTS/S/KeV')
    plt.title(f'{OBSID} Weighted I Spectra')
    #plt.ylim(0.01,10)
    plt.yscale('log')
    plt.legend()
    fig1 = plt.gcf()
    plt.show()
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/WeightedISpectra.pdf')
    plt.close()


#def WeightedQUspectra(OBSID,SOURCENAME,SOURCERAD):
    binwidth = 5/14 #The channel energy width is 0.04
    #ISpectra = plt.hist(df['E_AVG'], bins = np.arange(2,8,binwidth) )
    plt.figure()
    for det in range(1,3+1): #plus 1 so that we get all 3 detectors

        ISpectra = np.histogram(df['E_AVG'].loc[df['detector'] == det], bins = np.arange(2,8,binwidth)  )
        CenterofISpectraBins = (ISpectra[1][:-1] + ISpectra[1][1:])/2
        Q = []
        Q_Wmom = []
        W_mom = []
        Qstd = []
        I = []
        binwidth_2 = binwidth/2
        #for i in np.arange(2,8-binwidth,binwidth):
        for i in CenterofISpectraBins:
            I_ = len(df['Q'].loc[df['E_AVG'].between(i-binwidth_2,i+binwidth_2)& (df['detector'] == det)])
            Q_ = df['Q'].loc[df['E_AVG'].between(i-binwidth,i+binwidth_2)& (df['detector'] == det)].sum()
            Q_std = df['Q'].loc[df['E_AVG'].between(i-binwidth,i+binwidth_2)& (df['detector'] == det)].std()
            Q_Wmom_ = df['Q*W_mom'].loc[df['E_AVG'].between(i-binwidth_2,i+binwidth_2)& (df['detector'] == det)].sum()
            Wmom_ = df['W_mom'].loc[df['E_AVG'].between(i-binwidth_2,i+binwidth_2)& (df['detector'] == det)].mean()
            Q = np.append(Q,Q_)
            Q_Wmom = np.append(Q_Wmom, Q_Wmom_)
            W_mom = np.append(W_mom,Wmom_)
            Qstd = np.append(Qstd,Q_std)
            I = np.append(I,I_)
        Counts_Sec_KEV = (Q/LIVETIME/binwidth)
        WeightedCounts_Sec_KEV = (Q_Wmom/LIVETIME/binwidth)
        Q_Wmom_std = np.sqrt(abs(Qstd**2)*I)*W_mom/LIVETIME/binwidth

        plt.errorbar(CenterofISpectraBins,WeightedCounts_Sec_KEV,yerr = Q_Wmom_std,xerr = binwidth/2,markersize=3,label=f'Detector: {det}',fmt='o')
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Counts/S/KeV')
    plt.title(f'{OBSID} Weighted Q Spectra')
    plt.ylim(-0.004,0.004)
    plt.hlines(0,2,8,color = 'grey')
    plt.xlim(2,8)
    #plt.yscale('log')
    plt.legend()
    fig1 = plt.gcf()
    plt.show()
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/WeightedQSpectra.pdf')
    plt.close()

    #ISpectra = plt.hist(df['E_AVG'], bins = np.arange(2,8,binwidth) )
    plt.figure()
    for det in range(1,3+1): #plus 1 so that we get all 3 detectors

        ISpectra = np.histogram(df['E_AVG'].loc[df['detector'] == det], bins = np.arange(2,8,binwidth)  )
        CenterofISpectraBins = (ISpectra[1][:-1] + ISpectra[1][1:])/2
        I = []
        U = []
        U_Wmom = []
        W_mom = []
        Ustd = []
        for i in np.arange(2,8-binwidth,binwidth):
            I_ = len(df['U'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)])
            U_ = df['U'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].sum()
            U_std = df['U'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].std()
            U_Wmom_ = df['U*W_mom'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].sum()
            Wmom_ = df['W_mom'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].mean()
            U = np.append(U,U_)
            U_Wmom = np.append(U_Wmom, U_Wmom_)
            W_mom = np.append(W_mom,Wmom_)
            Ustd = np.append(Ustd,U_std)
            I = np.append(I,I_)
        Counts_Sec_KEV = (U/LIVETIME/binwidth)
        WeightedCounts_Sec_KEV = (U_Wmom/LIVETIME/binwidth)
        U_Wmom_std = np.sqrt(abs(Ustd**2)*I)*W_mom/LIVETIME/binwidth # I because that is the number of terms in our sum
        plt.errorbar(CenterofISpectraBins,WeightedCounts_Sec_KEV,yerr = U_Wmom_std,xerr = binwidth/2,markersize=3,label=f'Detector: {det}',fmt='o')
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Counts/S/KeV')
    plt.title(f'{OBSID} Weighted U Spectra')
    plt.ylim(-0.00,0.008)
    #plt.yscale('log')
    plt.legend()
    plt.hlines(0,2,8,color = 'grey')
    plt.xlim(2,8)
    fig1 = plt.gcf()
    plt.show()
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/WeightedUSpectra.pdf')
    plt.close()

#def NormalizedQUPlot(OBSID,SOURCENAME,SOURCERAD):
    binwidth = 5/14 #The channel energy width is 0.04
    #ISpectra = plt.hist(df['E_AVG'], bins = np.arange(2,8,binwidth) )
    plt.figure()
    for det in range(1,3+1): #plus 1 so that we get all 3 detectors    
        I_Wmom = []
        U_Wmom = []
        Energy = []
        Ustd = []
        I = []
        W_mom = []
        for i in np.arange(2,8,binwidth):
            CenterofISpectraBins_ = i +(binwidth/2) 
            I_ = df['W_mom'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].sum()
            U_Wmom_ = df['U*W_mom'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].sum()
            U_std = df['U'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].std()
            Wmom_ = df['W_mom'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)].mean()
            U_Wmom = np.append(U_Wmom, U_Wmom_)
            I_Wmom = np.append(I_Wmom, I_)
            Energy = np.append(Energy,CenterofISpectraBins_)
            Ustd = np.append(Ustd,U_std)
            W_mom = np.append(W_mom,Wmom_)
            I = np.append(I,len(df['U'].loc[df['E_AVG'].between(i,i+binwidth)& (df['detector'] == det)]))
        U_Wmom_std = np.sqrt(abs(Ustd**2)*I)
        U_Wmom_std = abs(np.sqrt(((abs(Ustd**2)*I*W_mom*W_mom)/(U_Wmom**2)) + ((I_Wmom**(1/4))/(I_Wmom**2)))*(U_Wmom/I_Wmom))


        plt.errorbar(Energy,U_Wmom/I_Wmom,markersize=3,yerr = U_Wmom_std,xerr = binwidth/2,label=f'Detector: {det}',fmt='o')
    plt.legend()
    plt.title(f'{OBSID} Normalized, Weighted U')
    plt.ylabel('U/I')
    plt.xlabel('Energy [KeV]')
    plt.hlines(0,2,8,color = 'grey')
    plt.xlim(2,8)
    fig1 = plt.gcf()
    plt.show()
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/NormalizedWeightedU.pdf')
    plt.close()

    Ebinwidth = binwidth #The channel energy width is 0.04
    #ISpectra = plt.hist(df['E_AVG'], bins = np.arange(2,8,binwidth) )
    plt.figure()
    for det in range(1,3+1): #plus 1 so that we get all 3 detectors    
        I_Wmom = []
        Q_Wmom = []
        Energy = []
        Qstd = []
        I = []
        W_mom = []
        for i in np.arange(2,8,Ebinwidth):
            CenterofISpectraBins_ = i +(Ebinwidth/2) 
            I_ = df['W_mom'].loc[df['E_AVG'].between(i,i+Ebinwidth)& (df['detector'] == det)].sum()
            Q_Wmom_ = df['Q*W_mom'].loc[df['E_AVG'].between(i,i+Ebinwidth)& (df['detector'] == det)].sum()
            Q_std = df['Q'].loc[df['E_AVG'].between(i,i+Ebinwidth)& (df['detector'] == det)].std()
            Wmom_ = df['W_mom'].loc[df['E_AVG'].between(i,i+Ebinwidth)& (df['detector'] == det)].mean()
            Q_Wmom = np.append(Q_Wmom, Q_Wmom_)
            I_Wmom = np.append(I_Wmom, I_)
            Energy = np.append(Energy,CenterofISpectraBins_)
            Qstd = np.append(Qstd,Q_std)
            W_mom = np.append(W_mom,Wmom_)
            I = np.append(I,len(df['Q'].loc[df['E_AVG'].between(i,i+Ebinwidth)& (df['detector'] == det)]))
        Q_Wmom_std = np.sqrt(abs(Qstd**2)*I)
        Q_Wmom_std = abs(np.sqrt(((abs(Qstd**2)*I*W_mom*W_mom)/(Q_Wmom**2)) + ((I_Wmom**(1/4))/(I_Wmom**2)))*(Q_Wmom/I_Wmom))


        plt.errorbar(Energy,Q_Wmom/I_Wmom,markersize=3,yerr = Q_Wmom_std,xerr = Ebinwidth/2,label=f'Detector: {det}',fmt='o')
    plt.legend()
    plt.title(f'{OBSID} Normalized, Weighted Q')
    plt.hlines(0,2,8,color = 'grey')
    plt.xlim(2,8)
    plt.ylabel('Q/I')
    plt.xlabel('Energy [KeV]')
    fig1 = plt.gcf()
    fig1.savefig(f'/home/jacob/Desktop/{SOURCENAME}/{OBSID}/plots/NormalizedWeightedQ.pdf')
    plt.show()
    #plt.savefig(f'/home/jacob/Desktop/HerX-1/{OBSID}/plots/NormalizedWeightedQ.pdf')
    plt.close()


#def QUPLOTS(OBSID,SOURCENAME,SOURCERAD):
    PBinWidth = 1/15
    PhaseBins = np.arange(0,1+PBinWidth,PBinWidth)
    #print(PhaseBins)
    U = []
    Q = []
    Ustd = []
    Qstd = []
    W_mom = []
    for i in np.arange(0,len(PhaseBins)-1,1):
        I = df['W_mom'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].sum()
        Wmom_ = df['W_mom'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].mean()

        U_ = df['U*W_mom'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].sum()
        U = np.append(U,U_/I)
        Ustd_ = df['U'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].std()
        Ustd = np.append(Ustd, np.sqrt((((Ustd_/I)**2))*I))


        Q_ = df['Q*W_mom'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].sum()
        Q = np.append(Q,Q_/I)
        Qstd_ = df['Q'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].std()

        Qstd = np.append(Qstd, np.sqrt((((Qstd_/I)**2))*I))

    fig,ax = plt.subplots(2,1,figsize=(5,7),gridspec_kw={'height_ratios':[3,1.5]})

    #Draw circles at 10 and 20% polarization
    circle = plt.Circle((0,0),0.1,fill = False,color = 'red',label ='10% Polarization')
    ax[0].add_patch(circle) 
    circle = plt.Circle((0,0),0.2,fill = False,color = 'green',label ='20% Polarization')
    ax[0].add_patch(circle)

    t = np.sin(PhaseBins[0:-1])
    ax[0].errorbar(Q,U,yerr = Ustd, xerr= Qstd,fmt='o',markersize = 0)
    ax[0].scatter(Q,U,s = 20,zorder=2,c=t)
    ax[0].set_ylim(-0.3,0.3)
    ax[0].set_xlim(-0.3,0.3)
    ax[0].hlines(0,-1,1,color = 'grey')
    ax[0].vlines(0,-1,1,color='grey')

    ax[0].set_xlabel('Q/I')
    ax[0].set_ylabel('U/I')
    ax[0].set_title(f'Phase Resolved U/I vs Q/I')
    #plt.legend()
    #plt.savefig(f'/home/jacob/Desktop/HerX-1/{OBSID}/plots/PhaseResolveQU.pdf')






    CenterofBins = (PulseP[1][:-1] + PulseP[1][1:])/2
    PULSEP = np.array([CenterofBins,PulseP[0]]).T
    PULSEP = pd.DataFrame(PULSEP, columns = ['Phase','Count'])
    Count = np.empty(0)
    CP = np.empty(0)
    for i in range(len(PhaseBins)-1):
        CP_ = (PhaseBins[i]+ PhaseBins[i+1])/2
        CP = np.append(CP,CP_)
        Count = np.append(Count,PULSEP['Count'].loc[PULSEP['Phase'].between(CP_-0.01,CP_+0.01)].mean())

    #plt.figure()
    ax[1].plot(CenterofBins,PulseP[0]/LIVETIME,color='red')

    ax[1].scatter(CP,Count/LIVETIME, c = t,zorder=2)
    ax[1].errorbar(CP,Count/LIVETIME,xerr = PBinWidth/2,color='blue',fmt='o',zorder=1)

    ax[1].set_xlabel('Phase')
    ax[1].tick_params(left = False,labelleft = False)
    ax[1].set_ylabel('Arbitrary Units')
    #ax[1].set_title(f'Pulse Profile')
    #plt.ylim(0,0.06)

    #ax[0].set_aspect('equal')
    #ax[1].set_aspect('equal')
    plt.tight_layout(pad = 1)
    fig1 = plt.gcf()
    plt.show
    fig1.savefig(f'/home/jacob/Desktop/HerX-1/{OBSID}/plots/QUPlot_PulseP.pdf')
    #plt.close()

    EBinWidth = 5/14
    U = []
    Q = []
    Ustd = []
    Qstd = []
    W_mom = []
    for i in np.arange(2,8,EBinWidth):
        I = df['W_mom'].loc[df['E_AVG'].between(i,i+EBinWidth)].sum()
        Wmom_ = df['W_mom'].loc[df['E_AVG'].between(i,i+EBinWidth)].mean()

        U_ = df['U*W_mom'].loc[df['E_AVG'].between(i,i+EBinWidth)].sum()
        U = np.append(U,U_/I)
        Ustd_ = df['U'].loc[df['E_AVG'].between(i,i+EBinWidth)].std()
        Ustd = np.append(Ustd, np.sqrt((((Ustd_/I)**2))*I))


        Q_ = df['Q*W_mom'].loc[df['E_AVG'].between(i,i+EBinWidth)].sum()
        Q = np.append(Q,Q_/I)
        Qstd_ = df['Q'].loc[df['E_AVG'].between(i,i+EBinWidth)].std()

        Qstd = np.append(Qstd, np.sqrt((((Qstd_/I)**2))*I))

    fig,ax = plt.subplots(2,1,figsize=(6,7))

    #Draw circles at 10 and 20% polarization
    circle = plt.Circle((0,0),0.1,fill = False,color = 'red',label ='10% Polarization')
    ax[0].add_patch(circle) 
    circle = plt.Circle((0,0),0.2,fill = False,color = 'green',label ='20% Polarization')
    ax[0].add_patch(circle)

    t = np.arange(2,8,EBinWidth)
    ax[0].errorbar(Q,U,yerr = Ustd, xerr= Qstd,fmt='o',markersize = 0)
    ax[0].scatter(Q,U,s = 20,zorder=2,c=t)
    ax[0].set_ylim(-0.3,0.3)
    ax[0].set_xlim(-0.3,0.3)
    ax[0].hlines(0,-1,1,color = 'grey')
    ax[0].vlines(0,-1,1,color='grey')

    ax[0].set_xlabel('Q/I')
    ax[0].set_ylabel('U/I')
    ax[0].set_title(f'Energy Resolved U/I vs Q/I')

    mappable= ax[1].scatter(np.arange(2,8,EBinWidth),np.full(len(np.empty_like(np.arange(2,8,EBinWidth))),0),c=t)
    #ax[1].set_xlabel('Energy [KeV]')
    plt.tight_layout(pad=1)

    children = fig1.get_children()[1]
    #cset = ax[0].contourf(Q,U,t,)
    fig.colorbar(mappable,ax=ax[0],label='Energy [KeV]')
    ax[1].set_visible(False)
    ax[0].set_aspect('equal',adjustable='box')
    fig1 = plt.gcf()
    plt.show()
    #plt.legend()
    fig1.savefig(f'/home/jacob/Desktop/HerX-1/{OBSID}/plots/EnergyResolveQU.pdf')
    plt.close()


#def PhaseResolvedPolarization(OBSID,SOURCENAME,SOURCERAD):
    df['U*W_mom*ModF'] = df['U*W_mom'] / df['SPECRESP']
    df['Q*W_mom*ModF'] = df['Q*W_mom'] / df['SPECRESP']

    PA = []
    PD = []
    U_ = []
    Q_ = []
    Ustd = []
    Qstd = []
    Phase = []
    PDE= np.empty(0)
    PAE = np.empty(0)

    for i in range(len(PhaseBins)-1):
        I = df['W_mom'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].sum()
        Wmom_ = df['W_mom'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].mean()


        respstd = df['SPECRESP'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].std()
        resp = df['SPECRESP'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].mean()


        U = df['U*W_mom*ModF'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].sum()/I
        Q = df['Q*W_mom*ModF'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].sum()/I

        U_std = df['U'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].std()/I
        Q_std = df['Q'].loc[df['Phase'].between(PhaseBins[i],PhaseBins[i+1])].std()/I

        PDe = np.sqrt((((U*U_std)**2)+((Q*Q_std)**2))
                                /((U**2)+((Q**2))))*100

        PAe = np.sqrt(((Q*U_std)**2)+(((U*Q_std)**2))
                               /(((U**2)+(Q**2))**2))*0.5*180/np.pi

        PA = np.append(PA,np.degrees((0.5)*np.arctan2(U,Q)))
        PD = np.append(PD,np.sqrt((U**2)+(Q**2))*100)
        PDE = np.append(PDE,PDe*Wmom_)
        PAE = np.append(PAE,PAe*Wmom_)
        Phase = np.append(Phase,(PhaseBins[i]+PhaseBins[i+1])/2)
        PDE = np.sqrt(abs(PD))
        PAE = np.sqrt(abs(PA))
        U_ = np.append(U_,U)
        Q_ = np.append(Q_,Q)
        Ustd = np.append(Ustd,np.sqrt(abs(U))*Wmom_*resp)
        Qstd = np.append(Qstd,np.sqrt(abs(Q))*Wmom_*resp)

    #Plot everything twice
    Phase = np.append(Phase,Phase+1)
    PA = np.append(PA,PA)
    PD = np.append(PD,PD)
    PDE = np.append(PDE,PDE)
    PAE = np.append(PAE,PAE)
    U_ = np.append(U_,U_)
    Q_ = np.append(Q_,Q_)
    Ustd = np.append(Ustd,Ustd)
    Qstd = np.append(Qstd,Qstd)


    fig,ax = plt.subplots(5,1,figsize=(5,25))

    ax[0].plot(CenterofBins,PulseP[0]/LIVETIME,color='blue')
    ax[0].plot(CenterofBins+1,PulseP[0]/LIVETIME,color='blue')

    ax[2].errorbar(Phase,U_,xerr=PBinWidth/2,yerr = Ustd,fmt='o',markersize=5)
    ax[2].set_title('U/I')
    ax[2].set_ylim(-0.2,0.4)

    ax[1].errorbar(Phase,Q_,xerr=PBinWidth/2,fmt='o',yerr=Qstd,markersize=5)
    ax[1].set_title('Q/I')
    ax[1].set_ylim(-0.3,0.3)

    ax[3].errorbar(Phase,PD,xerr=PBinWidth/2,yerr=PDE,fmt='o',markersize=5)
    ax[3].set_title('Polarization Degree')
    ax[3].set_ylim(5,35)



    ax[4].errorbar(Phase,PA,xerr=PBinWidth/2,yerr=PAE,fmt='o',markersize=5)
    ax[4].set_title('Polarization Angle')
    ax[4].set_ylim(20,80)

    plt.tight_layout(pad=1)
    fig1=plt.gcf()
    plt.show()
    fig1.savefig(f'/home/jacob/Desktop/HerX-1/{OBSID}/plots/Polarization&PulseProfile.pdf')

#def FullAnalysis(OBSID,SOURCENAME,SOURCERAD):
    #ImportOrbitalCorrectionData(OBSID,SOURCENAME,SOURCERAD)
    #ImportModulationFactor(OBSID,SOURCENAME,SOURCERAD)
    #ImportDataFiles(OBSID,SOURCENAME,SOURCERAD)
    #SetupBackGroundandEnergyBounds(OBSID,SOURCENAME,SOURCERAD)
    #BinaryCorrection(OBSID,SOURCENAME,SOURCERAD)
    #PeriodCorrection(OBSID,SOURCENAME,SOURCERAD)
    #PulseProfile(OBSID,SOURCENAME,SOURCERAD)
    #LightCurve(OBSID,SOURCENAME,SOURCERAD)
    #ImportARFMRF(OBSID,SOURCENAME,SOURCERAD)
    #ChanneltoEnergyandWeights(OBSID,SOURCENAME,SOURCERAD)
    #ISpectra(OBSID,SOURCENAME,SOURCERAD)
    #WeightedQUspectra(OBSID,SOURCENAME,SOURCERAD)
    #NormalizedQUPlot(OBSID,SOURCENAME,SOURCERAD)
    #PhaseResolvedPolarization(OBSID,SOURCENAME,SOURCERAD)

