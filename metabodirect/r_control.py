# exploration.py

import os
import subprocess
from shutil import which
from platform import system
# from loguru import logger


# --------------------------------------------------
def write_r_script(rscript, outdir, metadata_file, groups, norm_method='max'):
    """Write R scripts based on the provided Rscripts templates"""

    # Fixing the paths obtained in Windows, so they can work well in R and in Python
    current_dir = os.getcwd().replace('\\', '/')
    metadata_file = metadata_file.replace('\\', '/')
    metabo_home = os.path.split(os.path.realpath(__file__))[0].replace('\\', '/')

    # Modify the specified R script template with the desired values
    r_in = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'R_scripts_templates', rscript)
    r_in = open(r_in)
    new_name = rscript.replace('_template_2_groups', '').replace('_template', '')
    r_file = os.path.join(outdir, new_name)
    r_fh = open(r_file, 'w')
    for line in r_in:
        r_fh.write(line.replace('%currentdir%', current_dir)
                   .replace('%metadata%', metadata_file)
                   .replace('%Metabo_HOME%', metabo_home)
                   .replace('%outdir%', os.path.split(outdir)[0])
                   .replace('%group1%', groups[0])
                   .replace('%group2%', groups[1] if len(groups) == 2 else 'NULL')
                   .replace('%norm_method%', norm_method)
                   )
    r_fh.close()

    return r_file


# --------------------------------------------------
def run_r(rscript):
    """Execute R scripts"""

    # Check if the script is being run in Windows, or Linux and Mac so the correct Rscript executable is chosen
    rscript_exe = 'Rscript.exe' if system() == 'Windows' else 'Rscript'

    # Execute the R script
    if which(rscript_exe) is not None:
        r_cmd = [rscript_exe, rscript]
        p = subprocess.Popen(
            r_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = p.communicate()
        # if out:
        #     print(out)

        # if err:
        #     logger.warning('An error has occurred while running R, '
        #                    'please see the error below')
        #     logger.error(err)
    else:
        print('Rscript not accesible')
