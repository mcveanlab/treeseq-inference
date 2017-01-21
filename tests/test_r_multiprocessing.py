import os
import multiprocessing
import time
num_processes = 3
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def test_r_process(pause):
    importr("pryr")
    n_called = robjects.r("times.called <- times.called+1")[0]
    if pause<1:
        #allocate some memory to the R subprocess, just to test
        robjects.r("test <- runif(10000000)")
        robjects.r("test <- test + pi")
    r_pid = robjects.r("Sys.getpid()")[0]
    r_mem = robjects.r("mem_used()")
    print("R process for worker {} is {} (Mem: {}). Pausing for {} seconds.".format(
        os.getpid(), r_pid, r_mem, pause))
    time.sleep(pause)
    return(r_pid, n_called)
    
    
pause_secs = [2,4,3,6,0.5,2,3,1,5,2,3,3]
results = {}
robjects.r("times.called <- 0")
with multiprocessing.Pool(processes=num_processes) as pool:
    for proc, n_called in pool.imap_unordered(test_r_process, pause_secs):
        results[proc]=max(n_called, results.get(proc) or 0)
print("The test function should have been called {} times".format(len(pause_secs)))
for pid,called in results.items():
    print("R process {} was called {} times".format(pid,called))