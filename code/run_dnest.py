import shutil
import subprocess
import time as tsys
import numpy as np
import copy
import glob
import argparse

import dnest4

def rewrite_main(filename, dnest_dir = "./"):
    '''Rewrite the main.cpp to include
    the correct path to the filename'''

    mfile = open(dnest_dir+"main.cpp", "r")
    mdata = mfile.readlines()
    mfile.close()

    ## replace filename in appropriate line:
    mdata[-6] = '\tData::get_instance().load("%s");\n'%filename

    mfile.close()

    mwrite_file = open(dnest_dir+"main.cpp.tmp", "w")

    for l in mdata:
        mwrite_file.write(l)

    mwrite_file.close()

    shutil.move(dnest_dir+"main.cpp.tmp", dnest_dir+"main.cpp")

    return

def rewrite_options(nlevels=300, dnest_dir="./"):

    mfile = open(dnest_dir+"OPTIONS", "r")
    mdata = mfile.readlines()
    mfile.close()
    
    mdata[-4] = '%i\t# maximum number of levels\n'%nlevels

    mwrite_file = open(dnest_dir+"OPTIONS.tmp", "w")

    for l in mdata:
        mwrite_file.write(l)

    mwrite_file.close()

    shutil.move(dnest_dir+"OPTIONS.tmp", dnest_dir+"OPTIONS")

    return


def remake_model(dnest_dir="./"):
    ''' Runs make '''
    tstart = tsys.process_time() # changed from clock to process_time
    subprocess.call(["make", "-C", dnest_dir])
    tsys.sleep(15)
    tend = tsys.process_time() # changed from clock to process_time

    return

# # Added this extra functions from: https://github.com/eggplantbren/DNest3/blob/master/postprocess.py
# # needed for postprocess_new
# def logdiffexp(x1, x2):
# 	biggest = x1
# 	xx1 = x1 - biggest
# 	xx2 = x2 - biggest
# 	result = np.log(np.exp(xx1) - np.exp(xx2)) + biggest
# 	return result

# # Added this extra functions from: https://github.com/eggplantbren/DNest3/blob/master/postprocess.py
# # needed for postprocess_new
# def logsumexp(values):
# 	biggest = np.max(values)
# 	x = values - biggest
# 	result = np.log(np.sum(np.exp(x))) + biggest
# 	return result

# def postprocess_new(temperature=1., numResampleLogX=1, plot=False, save_posterior=False):

#     cut = 0

#     try:
#         levels = np.atleast_2d(np.loadtxt("levels.txt"))
#         sample_info = np.atleast_2d(np.loadtxt("sample_info.txt"))
#         sample = np.atleast_2d(np.loadtxt("csample.txt"))
#     except IOError:
#         return None, None

#     sample = sample[int(cut*sample.shape[0]):, :]
#     sample_info = sample_info[int(cut*sample_info.shape[0]):, :]
   
#     if sample.shape[0] != sample_info.shape[0]:
#         print('# Size mismatch. Truncating...')
#         lowest = np.min([sample.shape[0], sample_info.shape[0]])
#         sample = sample[0:lowest, :]
#         sample_info = sample_info[0:lowest, :]

#     # Convert to lists of tuples
#     logl_levels = [(levels[i,1], levels[i, 2]) for i in range(0, levels.shape[0])] # logl, tiebreaker
#     logl_samples = [(sample_info[i, 1], sample_info[i, 2], i) for i in range(0, sample.shape[0])] # logl, tiebreaker, id
#     logx_samples = np.zeros((sample_info.shape[0], numResampleLogX))
#     logp_samples = np.zeros((sample_info.shape[0], numResampleLogX))
#     logP_samples = np.zeros((sample_info.shape[0], numResampleLogX))
#     P_samples = np.zeros((sample_info.shape[0], numResampleLogX))
#     logz_estimates = np.zeros((numResampleLogX, 1))
#     H_estimates = np.zeros((numResampleLogX, 1))

#     # Find sandwiching level for each sample
#     sandwich = sample_info[:,0].copy().astype('int')
#     for i in range(0, sample.shape[0]):
#         while sandwich[i] < levels.shape[0]-1 and logl_samples[i] > logl_levels[sandwich[i] + 1]:
#             sandwich[i] += 1

#     for z in range(0, numResampleLogX):
#         # For each level
#         for i in range(0, levels.shape[0]):
#             # Find the samples sandwiched by this level
#             which = np.nonzero(sandwich == i)[0]
#             logl_samples_thisLevel = [] # (logl, tieBreaker, ID)
#             for j in range(0, len(which)):
#                 logl_samples_thisLevel.append(copy.deepcopy(logl_samples[which[j]]))
#             logl_samples_thisLevel = sorted(logl_samples_thisLevel)
#             N = len(logl_samples_thisLevel)

#             # Generate intermediate logx values
#             logx_max = levels[i, 0]
#             if i == levels.shape[0]-1:
#                 logx_min = -1E300
#             else:
#                 logx_min = levels[i+1, 0]
#             Umin = np.exp(logx_min - logx_max)

#             if N == 0 or numResampleLogX > 1:
#                 U = Umin + (1. - Umin)*np.random.rand(len(which))
#             else:
#                 U = Umin + (1. - Umin)*np.linspace(1./(N+1), 1. - 1./(N+1), N)
#             logx_samples_thisLevel = np.sort(logx_max + np.log(U))[::-1]

#             for j in range(0, which.size):
#                 logx_samples[logl_samples_thisLevel[j][2]][z] = logx_samples_thisLevel[j]

#                 if j != which.size - 1:
#                     left = logx_samples_thisLevel[j+1]
#                 elif i == levels.shape[0]-1:
#                     left = -1E300
#                 else:
#                     left = levels[i+1][0]

#                 if j != 0:
#                     right = logx_samples_thisLevel[j-1]
#                 else:
#                     right = levels[i][0]

#                 logp_samples[logl_samples_thisLevel[j][2]][z] = np.log(0.5) + logdiffexp(right, left)

#         logl = sample_info[:,1]/temperature

#         logp_samples[:,z] = logp_samples[:,z] - logsumexp(logp_samples[:,z])
#         logP_samples[:,z] = logp_samples[:,z] + logl
#         logz_estimates[z] = logsumexp(logP_samples[:,z])
#         logP_samples[:,z] -= logz_estimates[z]
#         P_samples[:,z] = np.exp(logP_samples[:,z])
#         H_estimates[z] = -logz_estimates[z] + np.sum(P_samples[:,z]*logl)

#     if save_posterior:

#         P_samples = np.mean(P_samples, 1)
#         P_samples = P_samples/np.sum(P_samples)
#         logz_estimate = np.mean(logz_estimates)
#         logz_error = np.std(logz_estimates)
#         H_estimate = np.mean(H_estimates)
#         H_error = np.std(H_estimates)
#         ESS = np.exp(-np.sum(P_samples*np.log(P_samples+1E-300)))

#         print("log(Z) = " + str(logz_estimate) + " +- " + str(logz_error))
#         print("Information = " + str(H_estimate) + " +- " + str(H_error) + " nats.")
#         print("Effective sample size = " + str(ESS))

#         # Resample to uniform weight
#         N = int(ESS)
#         posterior_sample = np.zeros((N, sample.shape[1]))
#         w = P_samples
#         w = w/np.max(w)
#         np.savetxt('weights.txt', w) # Save weights
#         for i in range(0, N):
#             while True:
#                 which = np.random.randint(sample.shape[0])
#                 if np.random.rand() <= w[which]:
#                     break
#             posterior_sample[i,:] = sample[which,:]
#         np.savetxt("posterior_sample.txt", posterior_sample)

#     return logx_samples, P_samples

def find_weights(p_samples):

    print("max(p_samples): %f" %np.max(p_samples[-10:]))

    ### NOTE: logx_samples runs from 0 to -120, but I'm interested in the values of p_samples near the
    ### smallest values of X, so I need to look at the end of the list
    if np.max(p_samples[-10:]) < 1.0e-5:
        print("Returning True")
        return True
    else:
        print("Returning False")
        return False

def run_burst(filename, dnest_dir = "./", levelfilename=None, nsims=100):
    '''' Automatically runs the DNest4 for the given filename without having
    to put the correct commands in the command line. Also finds the correct 
    amount of levels by running DNest4 twice. 

    Parameters:
    filename: the entire path to the file
    nsims: minimum amount of samples you want
    '''

    ### first run: set levels to 200
    # automatically runs the functions above
    print("Rewriting DNest run file (main.cpp)")
    rewrite_main(filename, dnest_dir)
    rewrite_options(nlevels=300, dnest_dir=dnest_dir) 
    remake_model(dnest_dir)

    # splits up the path to the file into the directory and the filename
    fdir = filename.split("/")
    fname = fdir[-1]
    fdir = filename[:-len(fname)]
    print("directory: %s" %fdir)
    print("filename: %s" %fname)

    fsplit = fname.split("_")
    # froot = "%s%s_%s" %(fdir, fsplit[0], fsplit[1]) # used to be: %s/%s_%s" 
    froot = "%s%s_%s" %(fdir, fsplit[0], fsplit[1]) # used to be: %s/%s_%s" 
    # froot = "%s%s_%s/%s_%s" %(fdir, fsplit[0], fsplit[1], fsplit[0], fsplit[1]) # used to be: %s/%s_%s" + additional fsplit[0] and fsplit[1] to create extra folder p/burst analysed
    print("froot: " + str(froot))

    print("First run of DNest: Find number of levels")
    ## run DNest
    # run a program with modified scheduling priority of -19 (using nice).
    # Niceness values range from -20 (most favorable to the  process)  to  19  (least  favorable  to  the process).
    dnest_process = subprocess.Popen(["nice", "-19", "./main", "-t", "8"])

    # endflag marks when DNest4 has to stop running 
    endflag = False
    while endflag is False:
        try:
            # every minute the plots pop up showing the log(x), iterations, etc. 
            tsys.sleep(60)
            logx_samples, p_samples = dnest4.postprocess(plot=False) # used to be postprocess_new
            # logx_samples, p_samples = postprocess_new(save_posterior=False)
            if p_samples is None:
                endflag = False
            else:
                endflag = find_weights(p_samples)
                print("Endflag: " + str(endflag))

        except (KeyboardInterrupt, ValueError):
            break
        
    # if the endflag is true print it 
    print("endflag: " + str(endflag))

    # Kill DNest4 when endflag is true: when weights are found 
    dnest_process.kill()

    ### CHECK FROM HERE ON: ##########
    dnest_data = np.loadtxt("%ssample.txt" %dnest_dir)
    nlevels = len(dnest_data)

    ### save levels to file
    if not levelfilename is None:
        levelfile = open(levelfilename, "a")
        levelfile.write("%s \t %i \n" %(filename, nlevels))
        levelfile.close()

    # Rerun DNest4 again, but now with the correct amount of levels calculated previously
    # For that rewrite the options file and the model(?)
    rewrite_options(nlevels=nlevels, dnest_dir=dnest_dir)
    remake_model(dnest_dir)

    dnest_process = subprocess.Popen(["nice", "-19", "./main", "-t", "8"])

    endflag = False
    while endflag is False:
        try:
            tsys.sleep(120)
            logx_samples, p_samples = dnest4.postprocess(plot=False) # used to be postprocess_new
            samples = np.loadtxt("%sposterior_sample.txt"%dnest_dir)
            print("samples file: %ssample.txt" %dnest_dir)
            print("nlevels: %i" %len(samples)) 
            print("Endflag: " + str(endflag))

            if len(samples) >= nsims and len(np.shape(samples)) > 1:
            #if len(samples) >= np.max([5*nlevels, 1000+nlevels]) and len(np.shape(samples)) > 1:
                endflag = True
            else:
                endflag = False
        except (KeyboardInterrupt, ValueError):
            break

    print("Endflag: " + str(endflag))

    # Kill DNest4 when endflag is true: this is now the case when ...?
    dnest_process.kill()

    # logx_samples, p_samples = dnest4.postprocess(save_posterior=True)  ##### DONT NEED IT??  
    # logx_samples, p_samples = postprocess_new(save_posterior=True)    

    # fsplit = filename.split("_") # DONE ALREADY BEFORE 
    # froot = "%s%s_%s" %(fdir, fsplit[0], fsplit[1][:-4]) ### used to be: "%s/%s_%s"
    # #froot = "%s_%s" %(fsplit[0], fsplit[1][-4])
    # print("froot: " + str(froot))

    # Move all the output to ...? 
    # Change the filenames so that it includes the filename(?)
    try:
        shutil.move("posterior_sample.txt", "%s_posterior_sample.txt" %froot)
        shutil.move("levels.txt", "%s_levels.txt" %froot)
        shutil.move("sample_info.txt", "%s_sample_info.txt" %froot)
        shutil.move("sample.txt", "%s_sample.txt" %froot)
        shutil.move("weights.txt", "%s_weights.txt" %froot)
        shutil.move("log_prior_weights.txt", "%s_log_prior_weights.txt" %froot) # added 
        shutil.move("main", "%s_main" %froot) # added 
        shutil.move("sampler_state.txt", "%s_sampler_state.txt" %froot) # added 

    except IOError:
        print("No file posterior_sample.txt")

    return

### NOT FINISHED YET ###
def run_all_bursts(data_dir="/home/mariska/UvA/magnetron2/data/", dnest_dir="./", levelfilename="test_levels.dat"):

    print("I am in run_all_bursts")
    # Run all the burst by running all the files that end with _test.dat
    filenames = glob.glob("%s*_test.dat"%data_dir)
    print(filenames)

    levelfilename = data_dir+levelfilename
    print("Saving levels in file %s"%levelfilename)

    levelfile = open(levelfilename, "w")
    levelfile.write("# data filename \t number of levels \n")
    levelfile.close()

    for f in filenames:
        print("Running on burst %s" %f)
        run_burst(f, dnest_dir=dnest_dir, levelfilename=levelfilename) 

    return

def main():
    print("I am in main")
    run_all_bursts(data_dir, dnest_dir, levelfilename)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Running DNest on a number of bursts")

    parser.add_argument("-d", "--datadir", action="store", required=False, dest="data_dir",
                        default="/home/mariska/UvA/magnetron2/data/", help="Specify directory with data files (default: current directory)")
    parser.add_argument("-n", "--dnestdir", action="store", required=False, dest="dnest_dir",
                        default="./", help="Specify directory with DNest model implementation "
                                           "(default: current directory")
    parser.add_argument("-f", "--filename", action="store", required=False, dest="filename",
                        default="test_levels.dat", help="Define filename for file that saves the number of levels to use")

    clargs = parser.parse_args()
    data_dir = clargs.data_dir
    dnest_dir = clargs.dnest_dir
    levelfilename = clargs.filename

    main()

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Running DNest on a number of bursts")

#     main(filename = "/home/mariska/UvA/magnetron2/data/B3_Jan14_test.dat")