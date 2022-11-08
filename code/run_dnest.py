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

def rewrite_options(nlevels=300, dnest_dir="./", output_path ="./"):

    mfile = open(dnest_dir+"OPTIONS", "r")
    mdata = mfile.readlines()
    mfile.close()
    
    mdata[-5] = '%i\t# maximum number of levels\n'%nlevels
    mdata[-1] = output_path

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

def run_burst(filename, dnest_dir = "./", fdir_out = "../output", levelfilename=None, nsims=100, min_nlevels=50): 
    '''' Automatically runs the DNest4 for the given filename without having
    to put the correct commands in the command line. Also finds the correct 
    amount of levels by running DNest4 twice. 

    Parameters:
    filename: the entire path to the file
    dnest_dir: relative path to the dnest directory 
    fdir_out: relative path to the folder with the output files 
    nsims: minimum amount of samples you want
    min_nlevels: minimum amount of levels yout want 
    '''

    # splits up the path to the file into the directory and the filename
    fdir = filename.split("/")
    fname = fdir[-1]
    fdir = filename[:-len(fname)]
    print("directory: %s" %fdir)
    print("filename: %s" %fname)

    fname = fname.removesuffix(".dat") 
    fsplit = fname.split("_")
        
    # make a seperate folder for the output of each filename 
    output_path = "%s/%s_%s/" %(fdir_out, fsplit[0], fsplit[1])
    froot= "%s/%s_%s/%s_%s" %(fdir_out, fsplit[0], fsplit[1], fsplit[0], fsplit[1]) 
    os.makedirs(os.path.dirname(output_path), exist_ok=True) 
    print("output_path: " + str(output_path))
    print("froot: " + str(froot))

    ### first DNest4 run: set levels to 300
    # automatically recompile and set the options file  
    print("Rewriting DNest run file (main.cpp)")
    rewrite_main(filename, dnest_dir)
    rewrite_options(nlevels=300, dnest_dir=dnest_dir, output_path=output_path) 
    remake_model(dnest_dir)

    print("First run of DNest: Find number of levels")
    # run a program with modified scheduling priority of -19 (using nice).
    # niceness values range from -20 (most favorable to the process) to 19 (least favorable to the process).
    dnest_process = subprocess.Popen(["nice", "-19", "./main", "-t", "8"])

    # endflag marks when DNest4 has to stop running 
    endflag = False
    while endflag is False:
        try:
            tsys.sleep(60)
            logz_estimate, H_estimate, logx_samples = dnest4.postprocess(plot=False, output_path=output_path)  
            levels_info = np.loadtxt("%slevels.txt" %output_path) 
            levels = len(levels_info)
            # set minimum amount of levels 
            if logx_samples is None or levels < min_nlevels : 
                endflag = False
            else:
                endflag = find_weights(logx_samples)
                print("Endflag: " + str(endflag))

        except (KeyboardInterrupt, ValueError):
            break
        
    # if the endflag is true print it 
    print("endflag: " + str(endflag))

    # Kill DNest4 when endflag is true: when weights are found 
    dnest_process.kill()

    dnest_data = np.loadtxt("%ssample.txt" %output_path)
    nlevels = len(dnest_data)

    ### save levels to file
    if not levelfilename is None:
        levelfile = open(levelfilename, "a")
        levelfile.write("%s \t %i \n" %(filename, nlevels))
        levelfile.close()

    # Rerun DNest4 again, but now with the correct amount of levels calculated previously
    # In order to do that rewrite the options file and recompile
    rewrite_options(nlevels=nlevels, dnest_dir=dnest_dir, output_path=output_path)
    remake_model(dnest_dir)

    dnest_process = subprocess.Popen(["nice", "-19", "./main", "-t", "8"])

    endflag = False
    while endflag is False:
        try:
            tsys.sleep(120)
            logz_estimate, H_estimate, logx_samples = dnest4.postprocess(plot=False, output_path=output_path) 
            samples = np.loadtxt("%sposterior_sample.txt"%output_path)
            print("samples file: %ssample.txt" %output_path)
            print("nsamples: %i" %len(samples)) 
            print("Endflag: " + str(endflag))

            if len(samples) >= nsims and len(np.shape(samples)) > 1: 
                endflag = True
            else:
                endflag = False

        except (KeyboardInterrupt, ValueError):
            break

    print("Endflag: " + str(endflag))

    # Kill DNest4 when endflag is true: 
    # when minimum amount of samples for the set amount of levels is reached 
    dnest_process.kill()

    # Change the outputfilenames so that it includes the filename
    try:
        shutil.move("%sposterior_sample.txt" %output_path, "%s_posterior_sample.txt" %froot) 
        shutil.move("%slevels.txt" %output_path, "%s_levels.txt" %froot)
        shutil.move("%ssample_info.txt" %output_path, "%s_sample_info.txt" %froot)
        shutil.move("%ssample.txt" %output_path, "%s_sample.txt" %froot)
        shutil.move("%sweights.txt" %output_path, "%s_weights.txt" %froot)
        shutil.move("%slog_prior_weights.txt" %output_path, "%s_log_prior_weights.txt" %froot)
        shutil.move("main", "%s_main" %froot)  # main is still not created in the output folder 

    except IOError:
        print("No file posterior_sample.txt")

    return

def run_all_bursts(data_dir="../data/", dnest_dir="./", fdir_out = "../output", levelfilename="test_levels.dat"):

    print("I am in run_all_bursts")
    # Run all the burst by running all the files that end with .dat
    filenames = glob.glob("%s*.dat"%data_dir) 
    print("Filenames:", filenames)

    levelfilename = data_dir+levelfilename
    print("Saving levels in file %s"%levelfilename)

    levelfile = open(levelfilename, "w")
    levelfile.write("# data filename \t number of levels \n")
    levelfile.close()

    for f in filenames:
        print("Running on burst %s" %f)
        run_burst(f, dnest_dir=dnest_dir, fdir_out=fdir_out, levelfilename=levelfilename) 

    return

def main():
    print("I am in main")
    run_all_bursts(data_dir, dnest_dir, fdir_out, levelfilename)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Running DNest on a number of bursts")

    parser.add_argument("-d", "--datadir", action="store", required=False, dest="data_dir",
                        default="../data/", help="Specify directory with data files (default: data directory)")
    parser.add_argument("-n", "--dnestdir", action="store", required=False, dest="dnest_dir",
                        default="./", help="Specify directory with DNest model implementation "
                                           "(default: current directory")
    parser.add_argument("--out", "--fdirout", action="store", required=False, dest="fdir_out",
                        default="../output", help="Specify directory for output files "
                                           "(default: output directory")
    parser.add_argument("-f", "--filename", action="store", required=False, dest="filename",
                        default="test_levels.dat", help="Define filename for file that saves the number of levels to use")

    clargs = parser.parse_args()
    data_dir = clargs.data_dir
    dnest_dir = clargs.dnest_dir
    fdir_out = clargs.fdir_out
    levelfilename = clargs.filename

    main()

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Running DNest on a single burst")

#     main(filename = "/home/mariska/UvA/magnetron2/data/B16_Jan14_profile2.dat")