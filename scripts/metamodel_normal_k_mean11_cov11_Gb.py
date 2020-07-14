import matlab.engine
import time
import sys
import os
import numpy as np

# set time
start=time.time()

# Gb_k
eng = matlab.engine.start_matlab()
eng.metamodel_normal_k_mean11_cov11_Gb(nargout=0)

eng.quit()

#running time
end=time.time()
print('running time={} s'.format(end-start))

