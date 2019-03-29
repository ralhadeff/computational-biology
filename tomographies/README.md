work in progress

the goal is to learn what cell slices look like and generate de-novo images of cell slices.

chronological order of attempts:  
`simple_AE-ul` - Autoencoder and generation from noise. Reconstruction seems good, but the image from noise gives the same averaged image.  
`vae_small` - VAE with no convolution, on a subset of the data.  
`manual_ae` - Keras AE and trying to input the hidden layer's input, using the hidden layer's distribution.  
`image_process` - more attention has been given to process the input images, removing noise and normalizing the intensity.  

So far all generated images are very blurry.

`greedy_layer` - greedy layer-wise pre-training. Normal AE is pre-trained on the data, and the hidden layer is used to run a VAE through. Some success, at least compared to previous attempts.
