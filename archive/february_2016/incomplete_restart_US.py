import os
import shutil

restart_path = "1_Restart"
os.makedirs(restart_path)

shutil.copy("run_npt.py",restart_path + "/run_npt.py")
shutil.copy("Ag_O1X5.5_O2X0.55.eam.fs",restart_path + "/Ag_O1X5.5_O2X0.55.eam.fs")
shutil.copy("pull_unix.in",restart_path + "/pull_unix.in")

