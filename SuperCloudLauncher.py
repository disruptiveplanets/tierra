import numpy as np
import os


BaseText = "#!/bin/bash \n#SBATCH -c 16 -n 1\n#SBATCH -o Run_%j.log\nsource /etc/profile\nmodule load anaconda/2022b\npython3 MCMC_Launch_WASP39.py"
#Instruments = ["NIRCAM","NIRISS","NIRSPEC_G395H","PRISM"]
Instruments = ["NIRSPEC_G395H"]

for CS in ["1", "4","5", "225", "2500"]:
    for Instr in Instruments:
        for PT_Profile in range(1,2): #(1,3) for parametric
            Command = BaseText+" "+CS+" "+Instr+" 50000 "+str(PT_Profile)
            with open("CurrentRun.sh", 'w') as f:
                f.write(Command)
            os.system("chmod u+x CurrentRun.sh")    
            os.system("LLsub CurrentRun.sh")
