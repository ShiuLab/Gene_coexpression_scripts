import os, sys

directory = sys.argv[1]
output_directory = "/mnt/home/uygunsah/4_Plastid/_June2016_resubmission/5_Bayesian/arc_lengths_random/done_so_far/" 
count=1
os.system('export PATH="/mnt/home/uygunsah/anaconda2/bin:$PATH"\n')
for file in os.listdir(directory):
    if file.startswith("selected"):
        #print file
	a = file.split("_")
	if len(a) == 5:
        	pathway1 = file.split("_")[3]
        	pathway2 = file.split("_")[4]
	else:
		pathway1 = file.split("_")[3]
		pathway2 = "100"
        naming = "arc_strength_pvalueMI_%s_%s" %(pathway1, pathway2)
        #if naming not in os.listdir(directory):
        #output =open("job%s.sh" %count, "w")
        #output.write("#!/bin/sh -login\n\n")
        #output.write("#PBS -q main\n")
        #output.write("#PBS -l nodes=1:ppn=1,walltime=3:59:00,mem=10gb\n")
        #output.write('export PATH="/mnt/home/uygunsah/anaconda2/bin:$PATH"\n')
        #output.write("R --vanilla --slave --args %s %s < bn_pathways.R" %(directory+file, naming))
        print "R --vanilla --slave --args %s %s < bn_pathways.R" %(directory+file, naming)
        os.system("R --vanilla --slave --args %s %s< bn_pathways.R" %(directory+file, naming))
        #output.close()
        
        #os.system("qsub job%s.sh" %count)
        #count = count+1
            
        
