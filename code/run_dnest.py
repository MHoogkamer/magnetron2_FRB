import shutil
import subprocess
import time as tsys
import numpy as np
import glob
import argparse
import os

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

    fname = fname.removesuffix(".dat") # ADDED 
    fsplit = fname.split("_")
        
    # make a seperate folder to put in the output for each filename 
    fdir_out = "/home/mariska/UvA/magnetron2/output_test"
    # froot = "%s%s_%s" %(fdir, fsplit[0], fsplit[1]) # used to be: %s/%s_%s" 
    froot= "%s/%s_%s/%s_%s" %(fdir_out, fsplit[0], fsplit[1], fsplit[0], fsplit[1])
    os.makedirs(os.path.dirname(froot), exist_ok=True)
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
            logx_samples, p_samples = dnest4.postprocess(plot=False) 
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

    dnest_data = np.loadtxt("%ssample.txt" %dnest_dir)
    nlevels = len(dnest_data)

    ### save levels to file
    if not levelfilename is None:
        levelfile = open(levelfilename, "a")
        levelfile.write("%s \t %i \n" %(filename, nlevels))
        levelfile.close()

    # Rerun DNest4 again, but now with the correct amount of levels calculated previously
    # In order to do that rewrite the options file and recompile
    rewrite_options(nlevels=nlevels, dnest_dir=dnest_dir)
    remake_model(dnest_dir)

    dnest_process = subprocess.Popen(["nice", "-19", "./main", "-t", "8"])

    endflag = False
    while endflag is False:
        try:
            tsys.sleep(120)
            logx_samples, p_samples = dnest4.postprocess(plot=False) 
            samples = np.loadtxt("%sposterior_sample.txt"%dnest_dir)
            print("samples file: %ssample.txt" %dnest_dir)
            print("nlevels: %i" %len(samples)) 
            print("Endflag: " + str(endflag))

            if len(samples) >= nsims and len(np.shape(samples)) > 1: 
                endflag = True
            else:
                endflag = False
        except (KeyboardInterrupt, ValueError):
            break

    print("Endflag: " + str(endflag))

    # Kill DNest4 when endflag is true: this is now the case when ...?
    dnest_process.kill()

    # Change the outputfilenames so that it includes the filename
    try:
        shutil.move("posterior_sample.txt", "%s_posterior_sample.txt" %froot)
        shutil.move("levels.txt", "%s_levels.txt" %froot)
        shutil.move("sample_info.txt", "%s_sample_info.txt" %froot)
        shutil.move("sample.txt", "%s_sample.txt" %froot)
        shutil.move("weights.txt", "%s_weights.txt" %froot)
        shutil.move("log_prior_weights.txt", "%s_log_prior_weights.txt" %froot) # added 
        shutil.move("main", "%s_main" %froot) # added 
        # shutil.move("sampler_state.txt", "%s_sampler_state.txt" %froot) # added 

    except IOError:
        print("No file posterior_sample.txt")

    return

def run_all_bursts(data_dir="/home/mariska/UvA/magnetron2/data/", dnest_dir="./", levelfilename="test_levels.dat"):

    print("I am in run_all_bursts")
    # Run all the burst by running all the files that end with .dat
    filenames = glob.glob("%s*.dat"%data_dir) ### NOTE: This has to be correct otherwise it doesn't work!
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