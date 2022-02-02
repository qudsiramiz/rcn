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
    #if len(file_list) >= 1501:
    #    # Check if the skip_rate is an integer
    #    if skip_rate != int(skip_rate):
    #        raise ValueError("skip_rate must be an integer")
    #    file_list = file_list[-1500::skip_rate]
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
                print(f"Processing image {count} of {len(file_list)}")
                img = iio.imread(filename)
                writer.append_data(img)
    writer.close()

    print(f"{vid_name} is created\n")


def make_gifs(
    img_folder="/home/cephadrius/Desktop/git/rcn/figures/all_ridge_plots/t96/gaussian_interpolation/",
    vid_folder="/home/cephadrius/Desktop/git/rcn/figures/moving_pictures/",
    vid_name="all_ridge_t96_gaussian",
    vid_type="mp4",
    skip_rate=1,
    duration=0.05,
    fps=25
    ):
    """
    Make gifs for the last n days. Default is 120 days, averaged over the last 30 days.

    Parameters
    ----------
    number_of_days : int, optional
        Number of days to be considered for plotting the gif. The default is 120.

    Returns
    -------
        None.
    """

    print(f"Code execution started at (UTC):" +
          f"{datetime.datetime.utcfromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}\n")

    file_list_dict = {}

    file_list_dict["file_list"] = np.sort(glob.glob(img_folder + "*.png"))

    for i,key in enumerate(list(file_list_dict.keys())):
        vid_name = f"{vid_folder}{vid_name}_{fps}fps.mp4"
        try:
            gif_maker(file_list_dict[key], vid_name, mode="I", skip_rate=skip_rate, vid_type=vid_type, fps=fps, duration=0.05)
        except ValueError as e:
            print(e)
            pass

image_inputs_t01 = {
    "img_folder":"/home/cephadrius/Desktop/git/rcn/figures/all_ridge_plots/t01/gaussian_interpolation/",
    "vid_folder":"/home/cephadrius/Desktop/git/rcn/figures/moving_pictures/",
    "vid_name":"all_ridge_t01_gaussian",
    "vid_type":"mp4",
    "skip_rate":1,
    "duration":0.05,
    "fps":1
}

image_inputs_t96 = {
    "img_folder":"/home/cephadrius/Desktop/git/rcn/figures/all_ridge_plots/t96/gaussian_interpolation/",
    "vid_folder":"/home/cephadrius/Desktop/git/rcn/figures/moving_pictures/",
    "vid_name":"all_ridge_t96_gaussian",
    "vid_type":"mp4",
    "skip_rate":1,
    "duration":0.05,
    "fps":1
}
make_gifs(**image_inputs_t96)