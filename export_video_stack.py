'''
Method for producing a video of a data collection 
    - not very robust 
    - requires manual specification of h5 files from VMXm eiger 
    - unsure if this would work with other beamlines?
'''

import h5py
import numpy as np
import imageio.v2 as imageio
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------
master_file = "/dls/mx/data/nr27313/nr27313-446/thaum/thaum_25_master.h5"       # master file
data_file   = "/dls/mx/data/nr27313/nr27313-446/thaum/thaum_25_000001.h5"       # data file with images
meta_file   = "/dls/mx/data/nr27313/nr27313-446/thaum/thaum_25_meta.h5"
output_video = "eiger_sequence.mp4"
brightness = 10
frames = ""
# ----------------------------------------------------

def plot_diagonal_intensity(frame, outfile=None):
    """
    Plot the intensity along the main diagonal of a 2D diffraction image.

    Parameters
    ----------
    frame : np.ndarray
        2D array representing a single diffraction image.
    outfile : str or None
        Optional filename to save the plot (e.g. 'diag_plot.png').
        If None, the plot is just shown.
    """

    # Ensure input is 2D
    if frame.ndim != 2:
        raise ValueError("Input frame must be a 2D array.")

    # Extract diagonal intensities
    diagonal = np.diag(frame)

    # Generate pixel index array
    pixels = np.arange(len(diagonal))

    # Plot
    plt.figure(figsize=(8, 5))
    plt.scatter(pixels, diagonal, lw=1.2)
    plt.title("Intensity Along Image Diagonal")
    plt.xlabel("Diagonal Pixel Index")
    plt.ylabel("Intensity (counts)")
    plt.grid(True, alpha=0.3)

    # Save or show
    if outfile:
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Diagonal intensity plot saved to: {outfile}")
    else:
        plt.show()


# ----------------------------------------------------
# 1. Read exposure time from the master file
# ----------------------------------------------------
with h5py.File(master_file, "r") as f_master:
    det = f_master["/entry/instrument/detector"]

    if "exposure_time" in det:
        exposure = float(det["exposure_time"][()])
    elif "count_time" in det:
        exposure = float(det["count_time"][()])
    else:
        raise RuntimeError("No exposure_time or count_time in master file")

fps = 1.0 / exposure
print(f"Exposure time: {exposure} s -> FPS: {fps}")


# ----------------------------------------------------
# 2. Load DECTRIS native corrections from meta file
# ----------------------------------------------------
with h5py.File(meta_file, "r") as f_meta:
    mask      = f_meta["mask"][()]        # 1 = bad pixel
    flatfield = f_meta["flatfield"][()]   # multiplicative
    offset    = f_meta["offset"][()]      # subtract this

print("Mask shape:", mask.shape)
print("Flatfield shape:", flatfield.shape)
print("Offset shape:", offset.shape)

# ----------------------------------------------------
# 3. Open the Eiger data file
# ----------------------------------------------------
with h5py.File(data_file, "r") as f_data:
    dset = f_data["/data"]            
    num_frames = dset.shape[0]
    print(f"Frames available: {num_frames}")
    print("Dataset shape:", dset.shape)

    if frames == "":
        frames = num_frames
    elif frames > num_frames + 1:
        print(f"Too many specified frames ({frames}); using max number in dataset ({num_frames})")
        frames = num_frames
    writer = imageio.get_writer(output_video, fps=fps, codec="libx264")

    for i in range(frames):
        print(f"Processing frame {i+1}/{frames}", end="\r")

        raw = dset[i].astype(np.float64)
        if i ==0 : plot_diagonal_intensity(raw, "raw")
        # ----------------------------------------------------
        # APPLY DETECTOR CORRECTIONS (matches DIALS + DECTRIS)
        # ----------------------------------------------------

        # 1) Subtract offset
        corrected = raw - offset[0]
        if i ==0 : plot_diagonal_intensity(corrected, "offset")

        # 2) Apply flatfield
        corrected *= flatfield
        if i ==0 : plot_diagonal_intensity(flatfield, "flatfield")

        # 3) Apply mask
        corrected[mask != 0] = 0.0
        if i ==0 : plot_diagonal_intensity(corrected, "mask")

        corrected = np.maximum(corrected, 1)
        log_img = np.log(corrected) * brightness
        if i ==0 : plot_diagonal_intensity(log_img, "log_image")

        # Auto contrast
        low, high = np.percentile(corrected, (0.5, 99.5))
        scaled = (corrected - low) / (high - low)
        if i ==0 : plot_diagonal_intensity(scaled, "scaled")

        scaled = np.clip(scaled, 0, 1)
        scaled = 1.0 - scaled  # invert
        if i ==0 : plot_diagonal_intensity(scaled, "inverted")

        frame8 = (scaled * 255).astype(np.uint8)
        writer.append_data(frame8)

    writer.close()

print("\nSaved video:", output_video)