""" Set of functions for converting microED data collected as compressed tiff files
    to MRC files for use in dials. These scripts can perform recentering as well as
    gain correction for data collected on the K2 without a beamstop.

    Marcus Gallagher-Jones
    marcus.gallagher-jones@rfi.ac.uk
    The Rosalind Franklin Institute
    2023/10/19"""

import tifffile
import numpy as np
import matplotlib.pyplot as plt
import mrcfile as mrc
from skimage.morphology import remove_small_holes
from skimage.morphology import remove_small_objects
import dm4
import argparse
import os
import glob
import traceback

def options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tiff_path", "--input_tiff_path", default="None", 
                        help="Path to the input tiff file")
    parser.add_argument("-defects", "--input_defect_file", default="None", 
                        help="Path to the input defects file")
    parser.add_argument("-gain_ref", "--input_gain_reference",default="None", 
                        help="Gain reference for subtraction from individual frames")
    parser.add_argument("-out", "--output_path", default="None", 
                        help="Output path for mrc file")
    parser.add_argument("-gain", "--gain_value", type=int, default=1, 
                        help="Gain value to apply to reduce rounding error")
    parser.add_argument("-binx", "--binning_factor_x", type=int, default=1, 
                        help="Amount to bin data in the x direction")
    parser.add_argument("-biny", "--binning_factor_y", type=int, default=1, 
                        help="Amount to bin data in the y direction")
    parser.add_argument("-binz", "--binning_factor_z", type=int, default=1, 
                        help="Amount to bin data in the z (time) direction")
    parser.add_argument("-thresh", "--centre_threshold", type=int, default=100, 
                        help="Intensity threshold for determening beam centre")
    parser.add_argument("-hot","--hot_pixel_val", type=int, default=1,
                        help="Threshold sigma level for removing hot pixels")
    parser.add_argument("-centre", "--perform_recentering", type=int, default=0, 
                        help="Binary flag to determine whether or not to recentre the diffraction images")
    parser.add_argument("-corr_gain", "--do_gain_correction", type=int, default=0, 
                        help="Binary flag to determine whether or not to gain correct the diffraction images")
    parser.add_argument("-batch", "--batch_size", type=int, default=0, 
                        help="Batch size to use when processing images with memmap")
    parser.add_argument("-pad", "--padded_size", type=int, default=0,
                        help="Size to pad the data up to such that recentering doesn't clip the data")
    parser.add_argument("-win", "--window_size", type=int, default=0,
                        help="Size of window for finding local hot pixels")
    parser.add_argument("-end", "--end_frame", type=int, default=[],
                        help="last frame to consider")
    args = parser.parse_args()
    return args

 
def create_circ_mask(yy, xx, imy, imx, radius):
    """Function for creating a circular mask embedded within a specific roi
       with a specified radius"""
    # Create coordinate grids for the y and x dimensions
    Y, X = np.meshgrid(np.arange(imx), np.arange(imy))
    
    # Calculate the squared distance from each point to the center of the circle
    distance_squared = (Y - xx)**2 + (X - yy)**2
    
    # Create a binary mask where pixels within the circle are set to 1, and others are 0
    circle_mask = (distance_squared <= radius**2)
    
    return circle_mask


def load_tif(filename, batch_size, batch_num):
    """Load LZW compressed tif file captured by serialEM."""
    print("Loading tif file ",filename, "\n")
    with tifffile.TiffFile(filename) as tif:
        batch_dims = np.arange(batch_num, batch_num+batch_size)
        print(batch_dims)
        data = tif.asarray(key=batch_dims)

    return data


def pad_data(datain, pady, padx):
    """Function to pad data with zeros by a specified ammount in each dimension"""
    # Only want to pad individual diffraction images so set z dim padding to 0
    print("Padding images\n")
    pad_widths = ((0, 0), (pady, pady), (padx, padx))
    padded_data = np.pad(datain, pad_widths, mode='constant', constant_values=0)

    return padded_data


def remove_defects(data_in, defects_file):
    """Load the defects .txt file and remove the bad pixels from each frame of the tif."""
    print("Removing defects\n")
    # First read the defects file line by line and extract the lines
    # that contain pixel coordinates. This will create a list of lists
    defect_coords_list = [line.strip().split(' ')[1:] 
                          for line in open(defects_file, 'r') 
                          if line[0] == 'B']
    # Once the lines have been seperated out extract flatten the list and convert
    # to a 2d array of coordinates
    flattened_list = [item for sublist in defect_coords_list for item in sublist]
    defect_coords_array = np.array(flattened_list,dtype=int).reshape(len(flattened_list)//2,2)
    # Divide coordinates by 2 for non super res data
    if data_in.shape[1] != 7676:
        defect_coords_array = defect_coords_array//2
    # Use the coordinates to set bad pixels to zero
    data_in[:,defect_coords_array[:,1],defect_coords_array[:,0]] = 0

    return data_in


def find_hot_pixels(data_in, thresh, threshold, window_size):
    """Function to find hot pixels on the detector outside those that are within
       the beam centre."""
    print("Removing hot pixels\n")
    # First take the sum of the input images to get a better idea of the location 
    # of hot pixels
    sum_data = np.sum(data_in,0)
    ydim, xdim = sum_data.shape
    hot_pix_mask = np.zeros((ydim, xdim))
    
    for ii in range(0,ydim,window_size):
        for jj in range(0,xdim,window_size):
            # Create a sliding window to assess the presence of hot pixels
            # in the image based on the threshold which represents the number
            # of standard deviations the pixel is above the mean of the local area
            
            subset = sum_data[ii:ii+window_size, jj:jj+window_size]
            std = np.std(subset)
            hot_pix_mask[ii:ii+window_size, jj:jj+window_size] = (subset > std * thresh).astype(int)
            
            
    
    #print(std)
    y0, x0 = find_beam_centre(sum_data, threshold, [])
    mask = create_circ_mask(y0, x0, ydim, xdim, 50)
    #hot_pix = sum_data > thresh*std
    hot_pix_mask[mask == 1] = 0
    # set hot pixels to zero
    hot_pix_inds = hot_pix_mask == 1
    data_in[:, hot_pix_inds] = 0

    return data_in
    
# Modified 2024/05/27 to account for the updated dm4files method of loading data

def load_gain(gain_file):
    
    # allow numpy gain files for downsampled or preprocessed gain references
    if gain_file.lower().endswith(".npy"):
        return np.load(gain_file)
        
    # First load in the gain data which will be in DM4 format
    with dm4.DM4File.open(gain_file) as gain_dm4:
        # Read the tags from the header to get the location of the image data and its dimensions
        tags = gain_dm4.read_directory()
        image_data_tag = tags.named_subdirs['ImageList'].unnamed_subdirs[1].named_subdirs['ImageData']
        image_tag = image_data_tag.named_tags['Data']
        XDim = gain_dm4.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[0])
        YDim = gain_dm4.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[1])

        # Read the data in as a numpy array and reshape it
        gain_image = np.array(gain_dm4.read_tag_data(image_tag), dtype=float)
        gain_image = np.reshape(gain_image, (YDim, XDim))

        # flip and rotate image so that it matches the rotation of the image data (this is from serialEM documentation)
        gain_image = np.rot90(np.fliplr(gain_image))
    
    return gain_image


def correct_gain(data_in, gain_image):
    """Load and apply the gain reference and rescale data """
    print("Performing Gain correction\n")
    zz, yy, xx = data_in.shape
    gain_corrected = np.zeros((zz, yy, xx))

    # Apply gain reference to all datasets in the stack
    for ii in range(0, zz):
        im = data_in[ii, :, :]
        im = im/gain_image
        im[im < 0] = 0
        gain_corrected[ii, :, :] = im
    
    return gain_corrected


def find_beam_centre(image, threshold, mask):
    """Function to find the centre of the beam using centre of mass"""
    if mask != []:
        image = image*mask
    std = np.std(image[image > 0])
    image_mask = image > int(threshold*std)
    image_mask = remove_small_holes(image_mask, 10000)
    image_mask = remove_small_objects(image_mask, 1000, connectivity=1)
     #plt.imshow(image_mask)
     #plt.show()
    image_mask = np.where(image_mask == 1)
    
    (y0, x0) = np.sum(image[image_mask] *image_mask, axis=1) / np.sum(image[image_mask])
    if np.isnan(y0) or np.isnan(x0):
        y0 = 0
        x0 = 0 
    return int(y0), int(x0)


def recentre_images(data_in, threshold, all_cents):
    """Function tpo  and then circshift image to centre"""
    print("Recentering images\n")
    zz, yy, xx = data_in.shape
    sum_data = np.sum(data_in,0)
     #y0, x0 = find_beam_centre(sum_data, threshold)
    for ii in range(0, zz):
        im = data_in[ii, :, :]
        y0, x0 = find_beam_centre(im, threshold, [])
        print(y0,x0)
        if y0 == 0:
            y0, x0 = all_cents[ii-1]
        all_cents.append((y0, x0))
        offsety = yy//2 - y0
        offsetx = xx//2 - x0
        im = np.roll(im, (offsety, offsetx), axis=(0, 1))
        data_in[ii, :, :] = im

    return data_in, all_cents


def bin_image_stack(data_in, binx, biny, binz):
    """Bin the raw frames by a set ammount in each dimension.
       Binning is done by summation to preserve intensity. Assumes
       that z is the first dimension"""
       
    # Reshape the data to create a higher dimensional array where every other 
    # axis is the local area over which to sum
    new_shape = [data_in.shape[0] // binz, 
                 data_in.shape[1] // binx, 
                 data_in.shape[2] // biny]

    # Check that rotation direction can actually be binned by the ammount specified
    # and crop the necessary frames if not
    difference = np.mod(data_in.shape[0],binz)
    if difference != 0:
        data_in = data_in[:-difference,:,:]

    print("Starting binning\n")
    new_layout = (new_shape[0], binz, new_shape[1], binx, new_shape[2], biny)
    print("Binning complete\n")
    print(new_shape[0])
    return data_in.reshape(new_layout).sum(-1).sum(-2).sum(1)


def write_2_mrc(data_in, output_path, scale_factor):
    """Convert the corrected tif stack to a mrc format that can be recognised by dials."""
    print("Writing mrc file ", output_path, "\n")
    data_in = data_in*scale_factor
    with mrc.new(output_path, overwrite=True) as out_file:
        out_file.set_data(np.uint16(data_in))


def convert_2_mrc(file_name, base_name, args):
    """Main function for converting the tif files into mrc forrmat"""
    
    # First open metadata of large tif file and figure out major dimensions
    tif_data = tifffile.TiffFile(file_name)
    nframes = len(tif_data.pages)
    print(f"nframes = {nframes}")
    frame_0 = tif_data.asarray(key=0)
    yy, xx = frame_0.shape
    pady = (args.padded_size - yy)//2
    padx = (args.padded_size - xx)//2
    if args.padded_size > 0:
        new_y = yy + (pady*2)
        new_x = xx + (padx*2)
    else:
        new_y = yy
        new_x = xx
    # Check if there are extra frames that do not fit neatly into the temporal
    # binning factor and create a holder for containing binned data
    extra_frames = nframes%args.binning_factor_z
    
    if args.end_frame:
        end_frame = args.end_frame
    else:
        end_frame = ((nframes-extra_frames)//args.binning_factor_z)
        
    resized_data = np.zeros((end_frame,
                        args.padded_size//args.binning_factor_y,
                        args.padded_size//args.binning_factor_x))
    if args.perform_recentering == 1:
        all_cents = []
    gain_image = load_gain(args.input_gain_reference)
    print("gain shape:", gain_image.shape)
    total_batch = (nframes-extra_frames)//args.batch_size
    print(f"total batch = {total_batch}")
    
   
    batch_ratio = args.batch_size//args.binning_factor_z
    # Loop through data in batches and create the final binned image    
    for ii in range(0, nframes - extra_frames, args.batch_size):
        batch_num = ii//args.batch_size
        if batch_num < total_batch:
            print(f"Analyzing batch number {batch_num} of {total_batch}\n")
            tif_data = load_tif(file_name, args.batch_size, ii)
            print(tif_data)
            tif_data = find_hot_pixels(tif_data, args.hot_pixel_val, args.centre_threshold, args.window_size)
            tif_data = remove_defects(tif_data, args.input_defect_file)
            if args.do_gain_correction == 1:
                tif_data = correct_gain(tif_data, gain_image)
            
            print("Frame min/max after gain:", 
                np.nanmin(tif_data[0]), np.nanmax(tif_data[0]))

            if args.padded_size > 0:
                tif_data = pad_data(tif_data, pady, padx)
            if args.perform_recentering == 1:
                tif_data, all_cents = recentre_images(tif_data, args.centre_threshold, all_cents)
            if batch_ratio > 1:
                resized_data[batch_num*batch_ratio:batch_num*batch_ratio + batch_ratio,:,:] = bin_image_stack(tif_data, 
                                                                                                            args.binning_factor_x, 
                                                                                                            args.binning_factor_y, 
                                                                                                            args.binning_factor_z)
            else:    
                if ii == 0:
                    resized_data[ii,:,:] = bin_image_stack(tif_data, 
                                                    args.binning_factor_x, 
                                                    args.binning_factor_y, 
                                                    args.binning_factor_z)
                else:
                    resized_data[ii//args.binning_factor_z,:,:] = bin_image_stack(tif_data, 
                                                                            args.binning_factor_x, 
                                                                            args.binning_factor_y, 
                                                                            args.binning_factor_z)
        else:
            print(f"Analyzing batch number {batch_num} of {total_batch}\n")
            #args.batch_size = args.binning_factor_z
            args.batch_size = nframes - ii
            print(nframes - ii)
            tif_data = load_tif(file_name, args.batch_size, ii)
            tif_data = find_hot_pixels(tif_data, args.hot_pixel_val, args.centre_threshold, args.window_size)
            tif_data = remove_defects(tif_data, args.input_defect_file)
            if args.do_gain_correction == 1:
                tif_data = correct_gain(tif_data, gain_image)
            tif_data = pad_data(tif_data, pady, padx)
            if args.perform_recentering == 1:
                tif_data, all_cents = recentre_images(tif_data, args.centre_threshold, all_cents)
            if batch_ratio > 1:
                resized_data[batch_num*batch_ratio:batch_num*batch_ratio + args.batch_size,:,:] = bin_image_stack(tif_data, 
                                                                                                            args.binning_factor_x, 
                                                                                                            args.binning_factor_y, 
                                                                                                            args.binning_factor_z)
            else:    
                if ii == 0:
                    resized_data[ii,:,:] = bin_image_stack(tif_data, 
                                                    args.binning_factor_x, 
                                                    args.binning_factor_y, 
                                                    args.binning_factor_z)
                else:
                    resized_data[ii//args.binning_factor_z,:,:] = bin_image_stack(tif_data, 
                                                                            args.binning_factor_x, 
                                                                            args.binning_factor_y, 
                                                                            args.binning_factor_z)
         
    write_2_mrc(resized_data, args.output_path + base_name + '.mrc', args.gain_value)
    
    if args.perform_recentering ==1:
        y, x = zip(*all_cents)
        plt.scatter(x,y, marker='.', color='r', s=10)
        
        #plt.xlim(0, 4000)
        #plt.ylim(0, 4000)
        plt.gca().set_aspect('equal')
        plt.savefig('centre_positions.png',dpi=300)


args = options()
base_dir = args.input_tiff_path
output_dir = args.output_path
tiffs = glob.glob(f'{base_dir}/*tif', recursive=True)
original_batch = args.batch_size
print(">>> STARTING SCRIPT", flush=True)
print(f"captured arguments: {args}")
print(f"base_dir: {base_dir}")
print(f"output_dir: {output_dir}")
print(f"tiffs: {tiffs}")
print(f"original_batch: {original_batch}")
for tiff_file in tiffs:
    path, name = os.path.split(tiff_file)
    base_name = name.split('_')[-2]
    print(f"Analysing dataset {base_name}")
    if not os.path.isfile(args.output_path + base_name + '.mrc'):
        args.batch_size = original_batch
        convert_2_mrc(tiff_file, base_name, args)
