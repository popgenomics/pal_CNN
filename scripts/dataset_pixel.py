import pandas as pd
import numpy as np
import tensorflow as tf
import h5py
import os
import argparse

def read_csv_dataset(out, pictures, color_channel):
    channel_to_str = {1:'grayscale', 3:'rgb', 4:'rgba'}
    x_train, y_train, x_val, y_val, x_test, y_test = [], [], [], [], [], []
    train_name, val_name, test_name = [], [], []
    for filename in os.listdir(pictures):
        scenario = filename.split('_')[0]
        scenario = int(scenario.split('S')[1]) - 1
        simulation = int(filename.split('_')[1])
        filepath = f'{pictures}/{filename}'
        print(filepath)
        print(simulation)
        image = tf.keras.utils.load_img(path=filepath, color_mode=channel_to_str[color_channel])
        if simulation <= 2499:
            train_name.append(filename)
            x_train.append(image)
            y_train.append(scenario)
        elif simulation > 2499 and simulation <=3499:
            val_name.append(filename)
            x_val.append(image)
            y_val.append(scenario)
        elif simulation > 3499:
            test_name.append(filename)
            x_test.append(image)
            y_test.append(scenario)

    with h5py.File(out, "w") as f:
        f.create_dataset("x_train", data=x_train)
        f.create_dataset("y_train", data=y_train)
        f.create_dataset("x_val", data=x_val)
        f.create_dataset("y_val", data=y_val)
        f.create_dataset("x_test", data=x_test)
        f.create_dataset("y_test", data=y_test)
        f.create_dataset('train_name', data=train_name)
        f.create_dataset('val_name', data=val_name)
        f.create_dataset('test_name', data=test_name)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', '-o', help='Output path and name', type=str, required=True, action='store')
    parser.add_argument('--images', '-i', help='The repertory with the pictures', type=str, required=True, action='store')
    parser.add_argument('--channel', '-c', help='Number of image\'s channels', type=int, required=True, choices={1, 3, 4}, action='store')
    args = parser.parse_args()
    read_csv_dataset(out=args.out, pictures=args.images, color_channel=args.channel)
