import os
import argparse
import bounding_box as bb


def main():
    """
    Read bounding boxes.
    """

    parser = argparse.ArgumentParser(description='Read bounding boxes.')
    parser.add_argument('input', type=str, help='Bounding box file.')

    args = parser.parse_args()
    if not os.path.exists(args.input):
        print('Input file does not exist.')
        exit(1)

    bounding_boxes = bb.read_bounding_boxes(args.input)
    for bounding_box in bounding_boxes:
        print(bounding_box)


if __name__ == '__main__':
    main()