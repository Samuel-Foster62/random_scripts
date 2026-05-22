from dxtbx.model.experiment_list import ExperimentListFactory
import numpy as np
from PIL import Image

# -----------------------------------
# USER SETTINGS
# -----------------------------------
#"/dls/mx/data/nr27313/nr27313-432/hins/25-3/1um/1um_l5_1.nxs"
#"/dls/mx/data/nr27313/nr27313-486/us17/us17_full_2.nxs" 
#"/dls/mx/data/nr27313/nr27313-487/lysozyme/gp2/x_1.nxs"

nexus_file = "/dls/mx/data/nr27313/nr27313-446/thaum/thaum_25.nxs" 
num_frames = 2
brightness = 10  # same as dials.image_viewer Brightness slider
# -----------------------------------

# Load Nexus/HDF5 experiment
expts = ExperimentListFactory.from_filenames([nexus_file])
imageset = expts.imagesets()[0]

# Grab first frame for shape
first = imageset.get_raw_data(0)[0].as_numpy_array()
accum = np.zeros_like(first, dtype=np.float64)

# Stack first N frames
for i in range(num_frames):
    img = imageset.get_raw_data(i)[0].as_numpy_array()
    accum += img

# DIALS viewer uses log scaling, but with brightness multiplier
# Avoid log(0)
accum = np.maximum(accum, 1)

# --------- DIALS-LIKE LOG SCALING ----------
log_img = np.log(accum)

# --------- APPLY BRIGHTNESS ----------
# This mimics the behaviour of the brightness slider:
# final_pixel = log(pixel) * brightness
log_img *= brightness

# --------- AUTO CONTRAST (0.5–99.5 percentile) ----------
low, high = np.percentile(log_img, (0.5, 99.5))
scaled = (log_img - low) / (high - low)
scaled = np.clip(scaled, 0, 1)

#---------- INVERT GRAYSCALE ----------
scaled = 1.0 - scaled

# --------- SAVE 16-bit PNG ----------
out = (scaled * 65535).astype(np.uint16)
Image.fromarray(out).save("single_image.png")

print("Saved: single_image.png")

