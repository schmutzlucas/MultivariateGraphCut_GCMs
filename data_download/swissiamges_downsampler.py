import os
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
import math


def downscale_image_realistic(image_path, target_resolution, original_resolution=0.1, sigma=1):
    """
    Downscale an image to simulate a lower-resolution sensor.

    Parameters:
    - image_path (str): Path to the input high-resolution image.
    - target_resolution (float): Desired resolution in meters (e.g., 0.3 for 30cm).
    - original_resolution (float): Original image resolution in meters (default 0.1 for 10cm).
    - sigma (float): Standard deviation for Gaussian kernel.

    Returns:
    - PIL.Image: Downscaled image.
    """
    # Open image
    img = Image.open(image_path).convert('RGB')
    img_np = np.array(img, dtype=np.float32)

    # Apply Gaussian blur to simulate the sensor's PSF
    blurred = gaussian_filter(img_np, sigma=(sigma, sigma, 0))

    # Calculate the scale factor
    scale_factor = target_resolution / original_resolution

    if scale_factor <= 1:
        raise ValueError("Target resolution must be higher (larger in meters) than original resolution.")

    # Calculate new dimensions
    new_width = math.floor(blurred.shape[1] / scale_factor)
    new_height = math.floor(blurred.shape[0] / scale_factor)

    # Resize the image using bicubic interpolation for smoothness
    downscaled_img = Image.fromarray(blurred.astype(np.uint8)).resize((new_width, new_height), Image.BICUBIC)

    return downscaled_img


def process_images_in_folder(input_folder, output_folder, resolutions, original_resolution=0.1):
    """
    Process all images in a folder and save the downscaled versions.

    Parameters:
    - input_folder (str): Path to the folder containing the high-resolution images.
    - output_folder (str): Path to the folder where the downscaled images will be saved.
    - resolutions (list of tuples): List of target resolutions in meters and corresponding sigma values.
    - original_resolution (float): Original image resolution in meters (default 0.1 for 10cm).
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Process each image in the folder
    for filename in os.listdir(input_folder):
        if filename.endswith(('.jpg', '.jpeg', '.png', '.tif')):
            image_path = os.path.join(input_folder, filename)
            for res, sigma in resolutions:
                downscaled_img = downscale_image_realistic(
                    image_path=image_path,
                    target_resolution=res,
                    original_resolution=original_resolution,
                    sigma=sigma
                )
                output_filename = f"{os.path.splitext(filename)[0]}_{int(res * 100)}cm_sigma{sigma}.tif"
                downscaled_img.save(os.path.join(output_folder, output_filename), format='TIFF')
                print(f"Saved {output_filename}")


# Example usage
if __name__ == "__main__":
    input_folder = 'swissimages_thun'
    output_folder = 'swissimages_thun_downscaled'

    # Define target resolutions in meters and corresponding sigma values
    target_resolutions = [
        (0.3, 1.0),  # 30cm with sigma ~1.0
        (0.5, 1.5),  # 50cm with sigma ~1.5
        (0.8, 1.9),
        (1.0, 2.5),  # 1m with sigma ~2.5
        (1.2, 3.0)  # 1.2m with sigma ~3.0
    ]

    # Process all images in the input folder
    process_images_in_folder(
        input_folder=input_folder,
        output_folder=output_folder,
        resolutions=target_resolutions,
        original_resolution=0.1  # 10cm
    )
