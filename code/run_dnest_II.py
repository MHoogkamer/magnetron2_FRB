import argparse
import glob
from multiprocessing import Pool, cpu_count
import os
from sys import exit
import shutil
import subprocess
import time as tsys
import numpy as np
import dnest4

def rewrite_options(nlevels=300, dnest_dir="./", output_path ="./", fname=""):

    mfile = open(dnest_dir+"OPTIONS", "r")
    mdata = mfile.readlines()
    mfile.close()
    
    mdata[-5] = '%i\t# maximum number of levels\n'%nlevels
    mdata[-1] = output_path

    mwrite_file = open(dnest_dir+fname+"_OPTIONS.tmp", "w")

    for l in mdata:
        mwrite_file.write(l)

    mwrite_file.close()

    shutil.move(dnest_dir+fname+"_OPTIONS.tmp", dnest_dir+fname+"_OPTIONS")

    return

def find_weights(p_samples):
    """logx_samples runs from 0 to -120, but I'm interested in the values of p_samples near the
    smallest values of X, so I need to look at the end of the list
    """
    p_sample_max = np.max(p_samples[-10:])
    print(f"Max p_sample value: {p_sample_max}")
    return p_sample_max < 1.0e-5

def run_burst(filename, dnest_dir, fdir_out, levelfilename=None, nsims=200, min_nlevels=50):  
    print("Running new burst")
    # splits up the path to the file into the directory and the filename
    path = filename.split("/")
    fname, fdir = path[-1], path[-2]

    print(f"> directory: ../{fdir}/")
    print(f"> filename: {fname}")

    fname = fname.removesuffix(".dat") 

    # make a seperate folder for the output of each filename 
    output_path = f"{fdir_out}/{fname}/"
    os.makedirs(output_path, exist_ok=True) 
    print(f"> output_path: {output_path}")

    rewrite_options(nlevels=500, dnest_dir=dnest_dir, output_path=output_path, fname=fname)  

    print("> starting 1st run DNEST4")
    print("> ./main", "-t", "8", filename)

    dnest_process = subprocess.Popen(["nice", "-19", "./main", "-t", "8", "-o", f"{fname}_OPTIONS" , filename]) 
    
    # marks when DNest4 has to stop running 
    while True:
        try:
            tsys.sleep(20)
            logz_estimate, H_estimate, logx_samples = dnest4.postprocess(plot=False, output_path=output_path)
            levels_info = np.loadtxt(f"{output_path}levels.txt") 
            levels = len(levels_info)

            if logx_samples is None or levels < min_nlevels:
                continue # keep searching

            found = find_weights(logx_samples)
            if found:
                print("Found weights!")
                break

        except (KeyboardInterrupt, ValueError):
            break

    dnest_process.kill() # stop dnest

    dnest_data = np.loadtxt(f"{output_path}sample.txt")
    nlevels = len(dnest_data)

    # save levels to file
    if not levelfilename is None:
        with open(levelfilename, 'a') as levelfile:
            levelfile.write(f"{filename}\t{nlevels}\n")

    # rerun DNest4 again, but now with the correct amount of levels calculated previously
    # in order to do that rewrite the options file 
    rewrite_options(nlevels=nlevels, dnest_dir=dnest_dir, output_path=output_path, fname=fname)

    dnest_process = subprocess.Popen(["nice", "-19", "./main", "-t", "8", "-o", f"{fname}_OPTIONS" , filename]) # make sure in main.cpp lin 11 you've: "Data::get_instance().load(argv[5]);"

    # marks when DNest4 has to stop running 
    while True:
        try:
            tsys.sleep(20)
            logz_estimate, H_estimate, logx_samples = dnest4.postprocess(plot=False, output_path=output_path)
            samples = np.loadtxt(f"{output_path}posterior_sample.txt")
            print(f"samples file: {output_path}sample.txt") 
            print(f"nsamples: {len(samples)}")

            # check if minimum amount of posterior samples is reached 
            if len(samples) >= nsims and len(np.shape(samples)) > 1:
                break # stop searching

        except (KeyboardInterrupt, ValueError):
            break
    
    # kill DNest4 when when minimum amount of samples for the set amount of levels is reached 
    dnest_process.kill()

    froot = f"{fdir_out}/{fname}/{fname}"

    try:
        # change the output filenames so that it includes the filename
        shutil.move(f"{output_path}posterior_sample.txt", f"{froot}_posterior_sample.txt") 
        # shutil.move(f"{output_path}levels.txt", f"{froot}_levels.txt")
        # shutil.move(f"{output_path}sample_info.txt", f"{froot}_sample_info.txt")
        # shutil.move(f"{output_path}sample.txt", f"{froot}_sample.txt")
        # shutil.move(f"{output_path}weights.txt", f"{froot}_weights.txt")
        # shutil.move(f"{output_path}log_prior_weights.txt", f"{froot}_log_prior_weights.txt")

        # remove unnecessary files  
        os.remove(f"{fname}_OPTIONS")
        os.remove(f"{output_path}levels.txt")
        os.remove(f"{output_path}sample_info.txt")
        os.remove(f"{output_path}sample.txt")
        os.remove(f"{output_path}weights.txt")
        os.remove(f"{output_path}log_prior_weights.txt")


    except IOError:
        print("No file posterior_sample.txt")


def run_all_bursts(data_dir, dnest_dir, fdir_out, levelfilename):
    levelfilename = f"{fdir_out}/{levelfilename}"
    print("Saving levels in file", levelfilename)

    if levelfilename is None: 
        with open(levelfilename, 'a') as levelfile:
            levelfile.write("# data filename \t number of levels \n")

    # Run all the burst by running all the files that end with .dat
    filenames = glob.glob(f"{data_dir}*.dat")
    print("Filenames:", filenames)

    # Create arguments per run burst process
    args_per_process = [(f, dnest_dir, fdir_out, levelfilename) for f in filenames]

    poolsize = cpu_count() - 1
    print(f"Starting run_burst pool of {poolsize}")
    with Pool(poolsize) as pool:
        pool.starmap(run_burst, args_per_process)
    print("Done running all bursts")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Running DNest on a number of bursts")

    parser.add_argument("-d", "--datadir", required=False, dest="data_dir",
                        default="../data/", help="Specify directory with data files (default: data directory)")
    parser.add_argument("-n", "--dnestdir", required=False, dest="dnest_dir",
                        default="./", help="Specify directory with DNest model implementation "
                                           "(default: current directory")
    parser.add_argument("--out", "--fdirout", required=False, dest="fdir_out",
                        default="../output", help="Specify directory for output files "
                                           "(default: output directory")
    parser.add_argument("-f", "--filename", required=False, dest="filename",
                        default="test_levels.dat", help="Define filename for file that saves the number of levels to use")

    clargs = parser.parse_args()
    data_dir = clargs.data_dir
    dnest_dir = clargs.dnest_dir
    fdir_out = os.path.abspath(clargs.fdir_out)
    levelfilename = clargs.filename

    if not os.path.isdir(data_dir):
        print(f"'{data_dir}' is not a directory. Please provide a valid directory with data files!")
        exit(1) # exit with error
    
    os.makedirs(fdir_out, exist_ok=True)
    
    run_all_bursts(data_dir, dnest_dir, fdir_out, levelfilename)
