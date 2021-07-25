import datetime
import os
import subprocess
from EFTGenReader.GenReader.make_html import make_html
from EFTGenReader.GenReader.utils import clean_dir,get_files,move_files

def make_plots(sample_loc):
    if not os.path.exists(sample_loc):
        print("This file does not exist, exiting:",sample_loc)
        return

    macro_str = "exploreTH1EFT.C(\"{infs}\")".format(infs=sample_loc)
    cmd = ['root','-b','-l','-q',macro_str]

    print 'Root Command: {0}'.format(' '.join(cmd))
    subprocess.check_call(['root','-b','-l','-q',macro_str])
    print ""


def main():

    file_loc = "/afs/crc.nd.edu/user/k/kmohrman/EFT/root_files/EFTMaodHists_output_tree.root"
    #file_loc = "output/EFTMaodHists_output_tree.root"

    make_plots(file_loc)

    # Make the outpout directory
    timestamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')
    output_dir = "/afs/crc.nd.edu/user/k/kmohrman/www/EFT/TopCoffea/testing/{time}_exploreTH1EFT_debug".format(time=timestamp_tag)
    #output_dir = "/afs/crc.nd.edu/user/k/kmohrman/www/EFT/TopCoffea/testing/histeft_rwgt_debugging/TH1EFT/rwgt_debug_{time}_sm".format(time=timestamp_tag)

    if output_dir == "":
        print("Please specify an outpout directory, exiting...")
        raise Exception
    else:
        if not os.path.exists(output_dir):
            os.mkdir(output_dir);

    imgs = get_files(".",targets=["^.*\.png$"])
    move_files(files=imgs,target=output_dir)
    make_html(output_dir)

if __name__ == "__main__":
    main()

