import rawpy
import cv2
import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, label
import matplotlib.pyplot as plt
import os
from skimage.io import imsave
from skimage.draw import disk
import time
import imageio


# -----------------------------------------------------------------------------
# Raw2 to white balanced array functions
# -----------------------------------------------------------------------------

def read_rw2_file(file_path):

    with rawpy.imread(file_path) as raw:
            rgb_image = raw.postprocess()
            
    rotated_image = np.rot90(rgb_image, 2)
        
    return rotated_image

def white_balance(image_data):

    grey_card_region = image_data[:400, :, :]
    mean_rgb = grey_card_region.mean(axis=(0, 1))
    color_scaling_factors = mean_rgb / mean_rgb.mean()
    white_balanced_image = image_data / color_scaling_factors

    grey_card_region_adjusted = white_balanced_image[-400:, :, :]
    current_brightness = grey_card_region_adjusted.mean()
    target_brightness = 118
    brightness_scaling_factor = target_brightness / current_brightness

    final_image = np.clip(white_balanced_image * brightness_scaling_factor, 0, 255).astype(np.uint8)

    return final_image

# -----------------------------------------------------------------------------
# Display functions
# -----------------------------------------------------------------------------

def display_image(image_data):
    plt.figure(dpi=1200)
    plt.imshow(image_data, interpolation='nearest')
    plt.axis('off')
    plt.show()
    
def display_grey_image(image_data):
    plt.figure(dpi=1200)
    plt.imshow(image_data, cmap='gray', interpolation='nearest')
    plt.axis('off')
    plt.show()
    
def display_image_channels(image_data):

    r_channel = image_data[:, :, 0]
    g_channel = image_data[:, :, 1]
    b_channel = image_data[:, :, 2]
    
    plt.figure(dpi=1200, figsize=(10, 4))
    
    plt.subplot(1, 3, 1)
    plt.imshow(r_channel, cmap='gray', interpolation='nearest')
    plt.axis('off')
    plt.title('Red Channel')
    
    plt.subplot(1, 3, 2)
    plt.imshow(g_channel, cmap='gray', interpolation='nearest')
    plt.axis('off')
    plt.title('Green Channel')
    
    plt.subplot(1, 3, 3)
    plt.imshow(b_channel, cmap='gray', interpolation='nearest')
    plt.axis('off')
    plt.title('Blue Channel')
    
    plt.tight_layout()
    plt.show()
    
def display_cropped_images(cropped_images):
    num_images = len(cropped_images)
    
    if num_images == 0:
        print("No images to display.")
        return

    cols = 4
    rows = (num_images // cols) + (num_images % cols > 0)

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3))

    axes = axes.flatten() if num_images > 1 else [axes]

    for ax, (name, img) in zip(axes, cropped_images.items()):
        ax.imshow(img, interpolation='nearest')
        ax.set_title(name, fontsize=10)
        ax.axis('off')

    for ax in axes[num_images:]:
        ax.axis('off')

    plt.tight_layout()
    plt.show()

    
# -----------------------------------------------------------------------------
# Pixel segmentation functions
# -----------------------------------------------------------------------------

def get_green_dominant_pixels(image_data):

    r_channel = image_data[:, :, 0]
    g_channel = image_data[:, :, 1]
    b_channel = image_data[:, :, 2]
    
    green_dominant = (g_channel > r_channel) & (g_channel > b_channel)
    
    return(green_dominant)
    
def get_galleria_pixels(image_data):

    r_channel = image_data[:, :, 0]
    g_channel = image_data[:, :, 1]
    b_channel = image_data[:, :, 2]
    
    red_dominant_or_low_green = ((r_channel > g_channel) & (r_channel > b_channel)) | (g_channel < 70)
    
    return red_dominant_or_low_green

# -----------------------------------------------------------------------------
# Morphological Functions
# -----------------------------------------------------------------------------


def hole_fill(array, iterations=10, n = 2):

    structure = np.ones((n,) * array.ndim, dtype=bool)

    expanded = binary_dilation(array, structure=structure, iterations=iterations)
    filled = binary_erosion(expanded, structure=structure, iterations=iterations)
    
    return filled

def hole_remove(array, structure=None, iterations=50, n = 2):

    if structure is None:
        structure = np.ones((n,) * array.ndim, dtype=bool)

    filled = binary_erosion(array, structure=structure, iterations=iterations)
    expanded = binary_dilation(filled, structure=structure, iterations=iterations)
    
    return expanded

def largest_white_object(binary_array):

    labeled_array, num_features = label(binary_array)
    
    if num_features == 0:
        return np.zeros_like(binary_array, dtype=bool)
    
    sizes = np.bincount(labeled_array.ravel())
    
    sizes[0] = 0
    
    largest_label = sizes.argmax()
    
    largest_object = (labeled_array == largest_label)
    
    return largest_object

def remove_small_holes(binary_array, size_threshold):

    inverted_array = np.logical_not(binary_array).astype(int)
    
    labeled_array, num_features = label(inverted_array)
    
    for label_num in range(1, num_features + 1):
        if np.sum(labeled_array == label_num) > size_threshold:
            inverted_array[labeled_array == label_num] = 0
    
    result_array = inverted_array + binary_array
    
    return result_array

# -----------------------------------------------------------------------------
# Find Circle Centres
# -----------------------------------------------------------------------------

def find_circle_center(image_mask, start_x, start_y):
    directions = [
        (0, -1), (0, 1), (-1, 0), (1, 0),
        (-1, -1), (-1, 1), (1, -1), (1, 1)
    ]
    
    x, y = start_x, start_y
    radius = 395

    for step in [20, 10, 5, 2, 1]:
        while True:
            circle_mask = np.zeros_like(image_mask, dtype=bool)
            rr, cc = disk((y, x), radius, shape=image_mask.shape)
            circle_mask[rr, cc] = True

            current_score = np.sum(image_mask & circle_mask)
            best_move = None
            best_score = current_score

            for dx, dy in directions:
                new_x, new_y = x + dx * step, y + dy * step

                if not (0 <= new_x < image_mask.shape[1] and 0 <= new_y < image_mask.shape[0]):
                    continue

                test_mask = np.zeros_like(image_mask, dtype=bool)
                rr, cc = disk((new_y, new_x), radius, shape=image_mask.shape)
                test_mask[rr, cc] = True

                new_score = np.sum(image_mask & test_mask)

                if new_score > best_score:
                    best_score = new_score
                    best_move = (new_x, new_y)

            if best_move:
                x, y = best_move
            else:
                break

    return x, y

def find_all_circles(image_mask, coordinates):
    found_centers = {}

    for name, (start_x, start_y) in coordinates.items():
        center_x, center_y = find_circle_center(image_mask, start_x, start_y)
        found_centers[name] = (center_x, center_y)

    return found_centers

# -----------------------------------------------------------------------------
# Mask Functions
# -----------------------------------------------------------------------------

def apply_circle_masks_and_crop(image_data, circles):
    height, width, _ = image_data.shape
    radius = 380

    cropped_images = {}
    y, x = np.ogrid[:height, :width]

    for name, (center_x, center_y) in circles.items():
        x_min = max(0, center_x - radius)
        x_max = min(width, center_x + radius)
        y_min = max(0, center_y - radius)
        y_max = min(height, center_y + radius)

        mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
        mask = np.repeat(mask[:, :, np.newaxis], 3, axis=2)

        green_area = np.zeros_like(image_data)
        green_area[:, :] = [0, 255, 0]

        masked_image = np.where(mask, image_data, green_area)

        cropped_image = masked_image[y_min:y_max, x_min:x_max]
        cropped_images[name] = cropped_image

    return cropped_images


# -----------------------------------------------------------------------------
# Validation Image Function
# -----------------------------------------------------------------------------

def generate_validation_image(original_image, centres, galleria_masks, valfile):
    grayscale_image = cv2.cvtColor(original_image, cv2.COLOR_RGB2GRAY)
    
    validation_image = np.stack([grayscale_image] * 3, axis=-1)

    ORANGE = np.array([243, 156, 18], dtype=np.uint8)
    BLUE = np.array([23, 255, 245], dtype=np.uint8)
    PURPLE = np.array([194, 24, 255], dtype=np.uint8)

    height, width = grayscale_image.shape
    y, x = np.ogrid[:height, :width]

    def add_overlay(image, mask, color, alpha=0.5):
        overlay = np.zeros_like(image)
        overlay[mask] = color
        return np.clip(image + overlay * alpha, 0, 255).astype(np.uint8)

    white_balance_mask = np.zeros((height, width), dtype=bool)
    white_balance_mask[:400, :] = True
    validation_image = add_overlay(validation_image, white_balance_mask, ORANGE, alpha=0.6)

    radius = 380
    final_circle_masks = {}

    for n, (cx, cy) in centres.items():
        circle_mask = (x - cx) ** 2 + (y - cy) ** 2 <= radius ** 2

        if n in galleria_masks:
            y_min, y_max = max(0, cy - radius), min(height, cy + radius)
            x_min, x_max = max(0, cx - radius), min(width, cx + radius)
            
            galleria_mask_resized = np.zeros((height, width), dtype=bool)
            galleria_mask_resized[y_min:y_max, x_min:x_max] = galleria_masks[n]
            
            circle_mask = circle_mask & ~galleria_mask_resized

        final_circle_masks[n] = circle_mask

    for mask in final_circle_masks.values():
        validation_image = add_overlay(validation_image, mask, BLUE, alpha=0.5)

    for n, (cx, cy) in centres.items():
        if n in galleria_masks:
            y_min, y_max = max(0, cy - radius), min(height, cy + radius)
            x_min, x_max = max(0, cx - radius), min(width, cx + radius)

            galleria_mask_resized = np.zeros((height, width), dtype=bool)
            galleria_mask_resized[y_min:y_max, x_min:x_max] = galleria_masks[n]

            validation_image = add_overlay(validation_image, galleria_mask_resized, PURPLE, alpha=0.5)

    validation_output_path = f"{valfile}_validation.jpg"
    imageio.imwrite(validation_output_path, validation_image.astype('uint8'), format='jpeg')
    
# -----------------------------------------------------------------------------
# Batch Process Function
# -----------------------------------------------------------------------------

def process_rw2_files(input_dir, outdir, valdir):

    coordinates = {
        "A1": (950, 1100), "A2": (1850, 1100), "A3": (2750, 1100), "A4": (3650, 1100),
        "B1": (950, 2050), "B2": (1850, 2050), "B3": (2750, 2050), "B4": (3650, 2050),
        "C1": (950, 3000), "C2": (1850, 3000), "C3": (2750, 3000), "C4": (3650, 3000)
    }

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(valdir, exist_ok=True)
    
    files = [f for f in os.listdir(input_dir) if "RW2" in f]

    i = 0
    
    times = []

    for file in files:
        if file.lower().endswith(".rw2"):
            
            start_time = time.time()
            i += 1
            
            print(f"Processing file: {i} of {len(files)}")
            
            file_path = os.path.join(input_dir, file)
            base_name = os.path.splitext(file)[0]
            outfile = os.path.join(outdir, base_name)
            valfile = os.path.join(valdir, base_name)

            image_data = read_rw2_file(file_path)
            image_data = white_balance(image_data)
            green_dominant_pixels = get_green_dominant_pixels(image_data)
            green_dominant_pixels = hole_remove(green_dominant_pixels, n=2, iterations=10)
            green_dominant_pixels = hole_fill(green_dominant_pixels, n=2, iterations=50)
            centres = find_all_circles(green_dominant_pixels, coordinates)
            cropped = apply_circle_masks_and_crop(image_data, centres)

            galleria_masks = {}

            for n, img in cropped.items():
                galleria_pixels = get_galleria_pixels(img)
                galleria_pixels = largest_white_object(galleria_pixels)
                galleria_pixels = remove_small_holes(galleria_pixels, 2000)

                galleria_masks[n] = galleria_pixels

                processed_image = (img * galleria_pixels[:, :, np.newaxis])

                output_path = f"{outfile}_{n}.tiff"
                imsave(output_path, processed_image.astype('uint8'))

            generate_validation_image(image_data, centres, galleria_masks, valfile)
            
            tot_time = time.time() - start_time
            
            times.append(tot_time)
            
            mean_time = sum(times) / len(times)
            est_time = (mean_time * (len(files) - i)) / 60
            
            print(f"Done. (Execution time: {tot_time:.2f} seconds, Remaining: {est_time:.2f} minutes)")


indir = r'C:\Users\Ryan\Desktop\test'
outdir = r'C:\Users\Ryan\Desktop\test\individuals'
valdir = r'C:\Users\Ryan\Desktop\test\validations'

process_rw2_files(indir, outdir, valdir)

