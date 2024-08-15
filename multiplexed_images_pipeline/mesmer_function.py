import numpy as np
def Denoise_img(img,percentile):
    for i in range(0, img.shape[0]):
        img[i][img[i] > np.percentile(img[i][img[i] != 0], percentile)] = np.percentile(img[i][img[i] != 0], percentile)
        img[i][img[i] < np.percentile(img[i][img[i] != 0], 100 - percentile)] = 0
    return(img)

def Generate_mesmer_input(image, marker_list, in_nuc_marker_kwargs, in_mem_marker_kwargs):
    
    from skimage.exposure import rescale_intensity
    
    marker_dict = {j:i for i,j in enumerate(marker_list, start=0)}
    nuc_marker_kwargs = {}
    mem_marker_kwargs = {}

    for num, nuc_key in enumerate(in_nuc_marker_kwargs, start=0):
        markers = in_nuc_marker_kwargs[nuc_key]
        nuc_marker_kwargs[num] = [marker_dict[marker] for marker in markers]

    for num, mem_key in enumerate(in_mem_marker_kwargs, start=0):
        markers = in_mem_marker_kwargs[mem_key]
        mem_marker_kwargs[num] = [marker_dict[marker] for marker in markers]    


    out = np.zeros((1, image.shape[1], image.shape[2], 2), dtype='float32')
    rescaled_nuc = np.zeros((len(nuc_marker_kwargs), image.shape[1], image.shape[2]), dtype='float32')
    rescaled_mem = np.zeros((len(mem_marker_kwargs), image.shape[1], image.shape[2]), dtype='float32')
    rescaled_image = np.zeros(image.shape, dtype='float32')

    for channel in range(image.shape[0]):
        rescaled_image[channel] = rescale_intensity(image[channel], out_range=(0.0, 1.0))

    for key_num in nuc_marker_kwargs:
        markers = nuc_marker_kwargs[key_num]
        rescaled_nuc[key_num] = np.sum(rescaled_image[markers], axis=0)

    for key_num in mem_marker_kwargs:
        markers = mem_marker_kwargs[key_num]
        rescaled_mem[key_num] = np.sum(rescaled_image[markers], axis=0)

    if rescaled_nuc.shape[0] == 1:
         out[..., 0] = rescaled_nuc[0]
    else:
        for channel in range(rescaled_nuc.shape[0]):
            rescaled_nuc[channel] = rescale_intensity(rescaled_nuc[channel], out_range=(0.0, 1.0))
        out[..., 0] = np.sum(rescaled_nuc, axis=0)

    if rescaled_mem.shape[0] == 1:
         out[..., 1] = rescaled_mem[0]
    else:
        for channel in range(rescaled_mem.shape[0]):
            rescaled_mem[channel] = rescale_intensity(rescaled_mem[channel], out_range=(0.0, 1.0))
        out[..., 1] = np.sum(rescaled_mem, axis=0)
        
    return out

def Plot_mesmer_input(img_mesmer,title=''):
    rgb_images=create_rgb_image(img_mesmer,channel_colors=['blue', 'green'])
    # select index for displaying
    idx = 0
    
    if (title==''):
        title='Membrane channel'

    # plot the data
    fig, ax = plt.subplots(1,3,figsize=(15, 15))
    ax[0].imshow(img_mesmer[idx, ..., 0])
    ax[1].imshow(img_mesmer[idx, ..., 1])
    ax[2].imshow(rgb_images[idx, ...])
    ax[0].set_title('Nuclear channel')
    ax[1].set_title(title)
    ax[2].set_title('Overlay')
    
    ax[0].get_xaxis().set_visible(False)
    ax[0].get_yaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    ax[1].get_yaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    ax[2].get_yaxis().set_visible(False)
    
    plt.show()
    
    
def Plot_segmentation_result(rgb_img,segmentation_predictions,idx=0):
    overlay_img=make_outline_overlay(rgb_data=rgb_img,predictions=segmentation_predictions)
    
    # select index for displaying
    x = 0

    # plot the data
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    ax.imshow(overlay_img[idx, ...])

    ax.set_title('Whole Cell Segmentation-Predictions')
    
    plt.axis('off')
    plt.show()