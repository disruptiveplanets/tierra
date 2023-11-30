import glob
import matplotlib.pyplot as plt
import h5py
import numpy as np
import corner
import os


#Now add the formatter of matplotlib
import matplotlib as mpl
mpl.rc('font', family='sans-serif', size=15)
mpl.rc('font', serif='Helvetica Neue')
mpl.rc('font', serif='Skia')
mpl.rc('text', usetex='True')
mpl.rc('ytick',**{'major.pad':5, 'color':'black', 'major.size':11,'major.width':1.5, 'minor.size':5,'minor.width':0.75})
mpl.rc('xtick',**{'major.pad':5, 'color':'black',  'major.size':11,'major.width':1.5, 'minor.size':5,'minor.width':0.75})
#mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('axes',**{'linewidth':1.0,'edgecolor':'navy'})


def GetParameterNames(FileName):
    #Construct the location to the file
    
    if "SERVER" in FileName.upper():
        Location = "ServerMCMCParams/BestParam"+FileName.split("/")[1].replace(".h5",".txt")
    else:    
        Location = "MCMCParams/BestParam"+FileName.split("/")[1].replace(".h5",".txt")

    if os.path.exists(Location):
        FileContent = open(Location).readlines()
        ParameterNames = [Item.split(":")[0] for Item in FileContent[1:]]
        for counter, Parameter in enumerate(ParameterNames):
            if "_" in Parameter:
                ParameterNames[counter] = "LogMR$_{\mathrm{"+Parameter.split("_")[1]+"}}$"
        return ParameterNames    
    else:
        return []
    


LocationProgressData = "ServerProgressData/"

AllFiles = glob.glob(LocationProgressData+"*.h5")
print("All the files are given by:", AllFiles)


#For the specific being run bere...
#AllFiles = glob.glob("ServerProgressData/RobinsonT78N1*1.h5")


for filename in AllFiles:
    print("Currently working with ", filename)
    
    with h5py.File(filename, "r") as f:
        SaveName = filename.split("/")[1].replace(".h5","")
        LogProbability = np.array(f['mcmc']['log_prob'])
        Chain = np.array(f['mcmc']['chain'])
        MeanLogProb1D = np.abs(np.mean(LogProbability, axis=1))
        SelectIndex = MeanLogProb1D>1e-4
        TotalNumberPoints = np.sum(SelectIndex)

        print("The total number of points is:", TotalNumberPoints)


        if np.sum(SelectIndex)<5:
            print("Files have not been properly saved yet. \n\n\n")
            continue

     
        ParameterNames = GetParameterNames(filename)

        LogProbability = LogProbability[SelectIndex,:]
        MeanLogProb1D = MeanLogProb1D[SelectIndex]

        print("The length of the selected index is", np.sum(SelectIndex))
        
        #Take the walker with the lowestMeanProbability to find the
        LowestMeanWalkerIndex = np.argmin(np.abs(np.mean(LogProbability, axis=0)))
        SelectedWalker = LogProbability[:,LowestMeanWalkerIndex]

        LogProbSTD = np.nanstd(SelectedWalker[-50:]) 
        MeanLogProb = np.abs(np.nanmean(SelectedWalker[-50:]))



        #Find the BurnIn 
        BurnIn = 0
        for i in range(2, len(MeanLogProb1D)-50):
            CurrentMean = abs(np.mean(SelectedWalker[i:i+50]))
            if CurrentMean<(MeanLogProb+2.0*LogProbSTD):
                BurnIn = i
                break     

        if TotalNumberPoints-BurnIn>10000:
            print("Use the last 10000 points to construct the burnin")
            BurnIn=TotalNumberPoints-10000        
        elif TotalNumberPoints>10000 and BurnIn<5000:
            print("Resetting  burnin as half.")        
            BurnIn=TotalNumberPoints//2
        elif TotalNumberPoints>4000 and BurnIn<2000:
            print("Resetting  burnin as 2000.")        
            BurnIn=2000

        print("Over-riding burnin as half the total number of points")
        BurnIn=TotalNumberPoints//2
        
       


        StepSize = np.arange(len(MeanLogProb1D))+1    

        #Mark the outlier walker     
        plt.figure(figsize=(8,6))
        _, Z = np.shape(LogProbability)
        for i in range(Z):
            plt.plot(StepSize, np.abs(LogProbability[:,i]))
        plt.plot(StepSize, MeanLogProb1D, "ko", zorder=100) 
        plt.axhline(y=MeanLogProb, color="cyan", lw="2", linestyle=":")    
        plt.axvline(x=BurnIn+1, color="red", lw="2", linestyle=":")      
        plt.xlabel("Step-Size", fontsize=20)
        plt.ylabel("Log Probability", fontsize=20)
        plt.ylim(MeanLogProb-2*LogProbSTD, MeanLogProb+100)
        plt.yscale('log')
        plt.tight_layout()
        if not("Server" in LocationProgressData):
            plt.savefig("CornerPlotFigures/%s_LogProbability_BeforeRemovingOutlier.png" %SaveName)
        elif not("Masked" in LocationProgressData):
            plt.savefig("ServerCornerPlotFigures/%s_LogProbability_BeforeRemovingOutlier.png" %SaveName)    
        else:
            plt.savefig("ServerFiguresMasked/%s_LogProbability_BeforeRemovingOutlier.png" %SaveName)        
        plt.close('all')    

        
        #Remove the unrun chain
        Chain = Chain[SelectIndex,:,:]

        #Remove the burn ins
        AllChain = Chain[BurnIn:,:,:]
        LogProbability2D = np.abs(LogProbability[BurnIn:,:].T)
       
        M,N = np.shape(LogProbability2D)
        MinLogProbability = np.min(LogProbability2D)
        WalkerMean = np.median(LogProbability2D, axis=1)
        #print("Change this to 20.0")
        BadWalkers = WalkerMean>MinLogProbability+15.0 #Change this to 20 later
        #BadWalkers = WalkerMean>MinLogProbability+100.0

        #Now find bad walker based on the amplitude of the change 
        WalkerSTD = np.std(LogProbability2D, axis=1)
        TotalSTD = np.std(LogProbability2D)

        NormalizedSTD = WalkerSTD/TotalSTD

        BadWalkersSTD = NormalizedSTD<1e-4   
        BadWalkers = np.logical_or(BadWalkers, BadWalkersSTD)     
        
        if np.sum(BadWalkers)>Z-2:
            print("Quite a number of bad walkers for this run...skipping this for now")
            continue

        #Create 2D map of the log probability
        LogProbability2DBadWalker = np.copy(LogProbability2D)
        LogProbability2DBadWalker[BadWalkers,:] = np.nan


        plt.figure(figsize=(8,6))
        plt.subplot(211)
        Image1 = plt.imshow( np.log10(LogProbability2D), aspect="auto", interpolation="None")
        plt.subplot(212)
        Image2 = plt.imshow(np.log10(LogProbability2DBadWalker), aspect="auto")
        plt.colorbar(Image2)
        plt.xlabel("Stepsize")
        plt.ylabel("Walker")
        if not("Server" in LocationProgressData):
            plt.savefig("CornerPlotFigures/%s_LogProbability2D.png" %SaveName)
        elif not("Masked" in LocationProgressData):
            plt.savefig("ServerCornerPlotFigures/%s_LogProbability2D.png" %SaveName)    
        else:  
            plt.savefig("ServerFiguresMasked/%s_LogProbability2D.png" %SaveName)
        plt.close('all')   


        #Remove the bad walker
        Chain = AllChain[:,~BadWalkers, :]
        X,Y,Z = np.shape(Chain) 
        FlatChain = Chain.reshape(X*Y, Z)

        #Convert the mass from the earth mass into Jupiter radius
        FlatChain[:,0] /= 317.907 
        

        #Now look at the mean log probability value
        MeanLogProbability = np.abs(np.mean(LogProbability, axis=1))
        


        #Make plot for the corner plot
        plt.figure(figsize=(12,8))

        if len(ParameterNames)>1:
            print("Generating corner plot WITH labels")
            corner.corner(FlatChain, labels=ParameterNames, title_quantiles = [0.16,0.5,0.84], show_titles=True)
            
        else:
            print("Generating corner plot WITHOUT labels")
            corner.corner(FlatChain, title_quantiles = [0.16,0.5,0.84], show_titles=True)
        #plt.tight_layout()
        if not("Server" in LocationProgressData):
            plt.savefig("CornerPlotFigures/%s_CornerPlot.png" %SaveName)
        elif not("Masked" in LocationProgressData):
            plt.savefig("ServerCornerPlotFigures/%s_CornerPlot.png" %SaveName)
        else:    
            plt.savefig("ServerFiguresMasked/%s_CornerPlot.png" %SaveName)
        plt.close('all')
        

        if not("Server" in LocationProgressData):
            MCMCSaveName = "MCMCData/"+SaveName+".npy"
        else:  
            MCMCSaveName = "ServerMCMCData/"+SaveName+".npy"
            
        print("Saving under", MCMCSaveName," and the shape of flat chain is ", np.shape(FlatChain))

        np.save(MCMCSaveName, FlatChain)
        XSteps = np.arange(X)+1

        #Making the walker data
        plt.figure(figsize=(8, 18))
        
        
        for i in range(1, Z+1):
            plt.subplot(Z, 1, i)
            for WalkerCounter, WalkerStatus in enumerate(BadWalkers):
                if not(WalkerStatus):
                    ColorName="black"
                else:
                    ColorName = "red"    
                plt.plot(XSteps, AllChain[:,WalkerCounter,i-1], color=ColorName)
            if len(ParameterNames)>1:
                plt.ylabel(ParameterNames[i-1])
        
        #xlabel for the last panel     
        plt.tight_layout()   
        if not("Server" in LocationProgressData):
            plt.savefig("Figures/%s_AllWalkers.png" %SaveName)
        elif not("Masked" in LocationProgressData):
            plt.savefig("ServerCornerPlotFigures/%s_AllWalkers.png" %SaveName)
        else:
            plt.savefig("ServerFiguresMasked/%s_AllWalkers.png" %SaveName)    
        plt.close('all')

        plt.figure(figsize=(8, 18))
        for i in range(1, Z+1):
            plt.subplot(Z, 1, i)
            for WalkerCounter, WalkerStatus in enumerate(BadWalkers):
                if not(WalkerStatus):  
                    plt.plot(XSteps, AllChain[:,WalkerCounter,i-1])
            if len(ParameterNames)>1:
                plt.ylabel(ParameterNames[i-1])
        
        #xlabel for the last panel     
        plt.tight_layout()   
        if not("Server" in LocationProgressData):
            plt.savefig("Figures/%s_GoodWalkers.png" %SaveName)
        elif not("Masked" in LocationProgressData):
            plt.savefig("ServerCornerPlotFigures/%s_GoodWalkers.png" %SaveName)
        else:  
            plt.savefig("ServerFiguresMasked/%s_GoodWalkers.png" %SaveName)
        plt.close('all')

    
