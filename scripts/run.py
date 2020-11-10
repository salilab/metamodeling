import matlab.engine
import time
import sys
import os
import numpy as np

# set time
start=time.time()

# --------------------
# run the individual model
# --------------------
# postprandial
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.postprandial_normal_basal(nargout=0)

eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.postprandial_normal(nargout=0)

eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.postprandial_t2d(nargout=0)

# pancreas
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.pancreas(nargout=0)

# exocytosis
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.exocytosis(nargout=0)

# signaling
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.signaling(nargout=0)

# metabolism
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.metabolism(nargout=0)

# GLP1R
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.screening(nargout=0)

# --------------------
# run the meta-model
# --------------------
# 420 inferences
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.metamodel_normal(nargout=0)

# GLP1
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.metamodel_normal_GLP1(nargout=0)

eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.metamodel_t2d_GLP1(nargout=0)

# incretin
eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.metamodel_normal_incretin(nargout=0)

eng = matlab.engine.start_matlab()
eng.add_bnet(nargout=0)
eng.metamodel_t2d_incretin(nargout=0)

eng.quit()

