#!/home/cephadrius/pyspedas/bin/python
# -*- coding: utf-8 -*-

import datetime
import glob as glob
import time

import imageio as iio
import numpy as np


def gif_maker(file_list, vid_name, mode="I", skip_rate=10, vid_type="mp4", duration=0.05, fps=25):
    """
    Make a gif from a list of images.

    Parameters
    ----------
    file_list : list
        List of image files.
    vid_name : str
        Name of the gif file.
    mode : str, optional
        Mode of the gif. The default is "I".
    skip_rate : int, optional
        Skip rate of the gif. The default is 10.
    vid_type : str, optional
        Type of the video. The default is "mp4".
    duration : float, optional
        Duration for which each image is displayed in gif. The default is 0.05.
    fps : int, optional
        Frames per second for mp4 video. The default is 25.

    Raises
    ------
    ValueError
        If the skip_rate is not an integer.
    ValueError
        If the duration is not a float.
    ValueError
        If the file_list is empty.
    ValueError
        If vid_name is empty.

    Returns
    -------
    None.
    """
    if file_list is None:
        raise ValueError("file_list is None")
    if vid_name is None:
        raise ValueError("vid_name is None. Please provide the name of the gif/video")
    if len(file_list) == 0:
        raise ValueError("file_list is empty")

    if vid_type == "gif":
        if duration != float(duration):
            raise ValueError("duration must be a float")
    if vid_type == "mp4":
        if fps != int(fps):
            raise ValueError("Frame rate (fps) must be an integer")

    count = 0
    if vid_type == "gif":
        with iio.get_writer(vid_name, mode=mode, duration=duration) as writer:
            for file in file_list:
                count += 1
                print(f"Processing image {count} of {len(file_list)}")
                image = iio.imread(file)
                writer.append_data(image)
    elif vid_type == "mp4":
        with iio.get_writer(vid_name, mode=mode, fps=fps) as writer:
            for filename in file_list:
                count += 1
                # Print the image number being processed every 10 images, in green color
                if count % 10 == 0:
                    print(f"Processed images ==> \033[92m {count} \033[00m of"
                          f" \033[91m {len(file_list)} \033[00m")
                img = iio.imread(filename)
                writer.append_data(img)
    writer.close()

    print(f"{vid_name} is created\n")


def make_gifs(
    img_folder="../figures/all_ridge_plots/t96/None_interpolation_mms3/time_series/",
    vid_folder="/home/cephadrius/Desktop/git/rcn/figures/moving_pictures/",
    vid_name="all_ridge_2hr_1min_vid",
    vid_type="mp4",
    skip_rate=1,
    duration=0.05,
    fps=10
):
    """
    Make a gif from a list of images.

    Parameters
    ----------
    img_folder : str, optional
        Folder containing the images. The default is
        "../figures/all_ridge_plots/t96/None_interpolation_mms3/time_series/".
    vid_folder : str, optional
        Folder to save the gif. The default is
        "/home/cephadrius/Desktop/git/rcn/figures/moving_pictures/".
    vid_name : str, optional
        Name of the gif file. The default is "all_ridge_2hr_1min_vid".
    vid_type : str, optional
        Type of the video. The default is "mp4".
    skip_rate : int, optional
        Skip rate of the gif. The default is 1.
    duration : float, optional
        Duration for which each image is displayed in gif. The default is 0.05.
    fps : int, optional
        Frames per second for mp4 video. The default is 10.

    Returns
    -------
        None.
    """

    print("Code execution started at (UTC):" +
          f"{datetime.datetime.utcfromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}\n")

    file_list_dict = {}

    file_list_dict["file_list"] = np.sort(glob.glob(img_folder + "*.png"))

    for i, key in enumerate(list(file_list_dict.keys())):
        vid_name = f"{vid_folder}{vid_name}_{fps}fps.mp4"
        try:
            gif_maker(file_list_dict[key], vid_name, mode="I", skip_rate=skip_rate,
                      vid_type=vid_type, fps=fps, duration=0.05)
        except ValueError as e:
            print(e)
            pass


# Find all the folders in the current directory
folder_list = glob.glob("../figures/all_ridge_plots/t96/None_interpolation_mms3/time_series_20170901_120000/")

for folder in folder_list[:1]:
    image_inputs_t96 = {
        "img_folder": folder,
        "vid_folder": "/home/cephadrius/Desktop/git/rcn/figures/moving_pictures/",
        "vid_name": folder.split("/")[-2],
        "vid_type": "mp4",
        "skip_rate": 1,
        "duration": 0.05,
        "fps": 10
    }
    print(f"Processing folder: {folder}")
    make_gifs(**image_inputs_t96)
