# Fast Fourier Transform Image Compression
Sequential and parallel fast Fourier transform algorithms for image compression.

- [Fast Fourier Transform Image Compression](#fast-fourier-transform-image-compression)
  * [Introduction](#introduction)
  * [Examples](#examples)
    + [`images/grig.jpg` (512x1024)](#-images-grigjpg---512x1024-)
    + [`images/graphics.jpg` (512x512)](#-images-graphicsjpg---512x512-)
    + [`images/philip.jpg` (2048x2048)](#-images-philipjpg---2048x2048-)
  * [Implementation Details](#implementation-details)
    + [Serial 1D FFT](#serial-1d-fft)
    + [Parallel 1D FFT](#parallel-1d-fft)
    + [Serial 2D FFT](#serial-2d-fft)
    + [Parallel 2D FFT](#parallel-2d-fft)
    + [Compression](#compression)
  * [Compilation & Execution](#compilation---execution)
    + [Requirements](#requirements)
    + [Example Usage](#example-usage)
  * [Analysis](#analysis)
  * [Conclusion](#conclusion)
  * [References](#references)

## Introduction
<p align="justify"> A fast Fourier transform algorithm quickly computes the discrete Fourier transform (DFT) of a sequence, or its inverse (IDFT). Fourier analysis converts a signal from its original domain (often time or space) to a representation in the frequency domain and vice versa. The DFT is obtained by decomposing a sequence of values into components of different frequencies. </p>

<p align="justify"> This project focuses on the application of fast Fourier transforms for image compression. A Fourier transform allows deconstructs an image to its frequncy domain. In this domain, one can remove insignificant frequncies, essentially stripping data from the image. One can then reconstruct the image with less frequncies to obtain a "compressed" image. For the sake of simplicity, image width and height are constrained to be a power of two and all images are greyscaled before fourier transform computations. </p>

In this project, four types algorithms were implemented and their execution times were analyzed based on the number of threads:
1. Serial 2D FFT w/ Serial 1D FFT
2. Parallel 2D FFT w/ Serial 1D FFT
3. Serial 2D FFT w/ Parallel 1D FFT
4. Parallel 2D FFT w/ Parallel 1D FFT

Note: **The program does not actually output a compressed jpeg image.** It only shows what a compressed image would look like since the ouput is always simply written to a png.

## Examples
### `images/grig.jpg` (512x1024)
Input | 10% | 2% | 1%
:---:|:---:|:---:|:---:|
![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/grig.jpg) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/grig_0.1.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/grig_0.02.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/grig_0.01.png) |

Notice the loss of detail, particularly in the girl's eyebrows, and the man's earring disapearing at higher rates of compression.

### `images/graphics.jpg` (512x512)
Input | 10% | 5% | 1%
:----:|:---:|:---:|:---:|
![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/graphics.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/graphics_0.1.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/graphics_0.05.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/graphics_0.01.png) |

Compression at lower rates removes some noise from this image (which is a render from [my ray tracer](https://github.com/joshuapjacob/computer-graphics)). However, compression at higher rates adds more noise.

### `images/philip.jpg` (2048x2048)
Input | 1% | 0.2% | 0.1%
:---:|:---:|:---:|:---:|
![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/philip.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/philip_0.01.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/philip_0.002.png) | ![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/compressed/philip_0.001.png) |

Notice the loss of detail, particularly in the hair.
## Implementation Details

All algorithms modify the data in-place.

### Serial 1D FFT
This is the standard Cooley-Tukey Radix-2 decimination-in-frequncy algorithm. This is a divide-and-conquer algorithm that recursively divides the transform into two pieces of size N/2 at each step, and is therefore limited to power-of-two sizes.

### Parallel 1D FFT
This is an optimized Cooley-Tukey FFT decimation-in-time algorithm that uses bit-reversal and a table for the exponentiation. Unfortunatly, the parallelization of my implementation here is quite limited and only involves breaking down loops over data into blocks for each thread.

### Serial 2D FFT
This is the standard sequantial 2D FFT that invloves perfoming a 1D FFT on each row of the image data matrix and then on each column of the resuling matrix.

### Parallel 2D FFT
This is a 2D FFT that makes each thread perform a 1D FFT on each row of a block of rows of the image data matrix. Then, after syncronization, each thread performs a 1D FFT on each column of a of the resulting matrix. When used with the parallel 1D FFT, the 1D FFT computations only use 2 threads. For example, if the this algorithm is used with 8 threads, the 2D computations are only done using 4 threads and since each 2D computation thread requires doing 1D computation, each of the 4 threads will only spawn 2 threads, making the total number of threads to be 8.

### Compression
After a fourier transform is done, the frequncies are sorted and only a certain fraction of frequncies are kept, the rest of the data in the matrix is set to 0. Finding exactly which data to keep is the often very slow with large data matrices. I'm aware that such precision is usually not neccesary and a randoming sampling to approximated the distribution of the frequncies and then choosing a cut-off would have been much more faster than finding the exact cut-off.

## Compilation & Execution
### Requirements
- [stb](https://github.com/nothings/stb) (to read/write images)
- [tbb](https://github.com/oneapi-src/oneTBB) (high level abstract threading for `std::sort`)

```
$ g++ -O3 main.cpp -D_GLIBCXX_PARALLEL -fopenmp -ltbb -std=c++17 -o main.out
$ ./main.out <image filename> <fraction of data to keep>
```
The output images are stored as `output1.png`, `output2.png`, `output3.png`, and `output4.png`. All outputs are essentially the same apart from minuscule differences caused by numerial precision differences of the fast Fourier transform algorithms.

### Example Usage
```
$ ./main.out images/philip.png 0.05 # keeps only 5% of the data in philip.png
```

## Analysis

Analysis of excution time was only performed for the largest image (`images/philip.png` 2048x2048). All tests were run on a machine with an eight core processor. Execution time is not the complete program execution time but only the execution time of the fourier transform function and its inverse.
  
![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/plots/execution_time.png)

![](https://raw.githubusercontent.com/joshuapjacob/fast-fourier-transform-image-compression/main/images/plots/speed_up.png)

The fully sequential algorithm (Serial 2D FFT w/ Serial 1D FFT) is independant of the number of threads as expected. The fully parallel algorithm (Parallel 2D FFT w/ Parallel 1D FFT) offers the best speed-up with a high number of threads. However, parallelizing only the 2D FFT (Parallel 2D FFT w/ Serial 1D) seems to offer a better speed-up when the number of threads is below 12. The worst algorithm involves parallelizing only the 1D FFT (Serial 2D FFT w/ Parallel 1D FFT). It does not offer as significant of a speed-up when the number of threads is over 3. In fact, it seems to progressivly get worse when the number of threads is increased. This might be due to thread-spawning overhead becoming significant or I just have a bug in my code :confused: .

## Conclusion

The overall image compression of all the algorithms implemented is successful. We observe that more information is lost when we increasing the compression rate. The best algorithm is the fully parallel implementation (Parallel 2D FFT w/ Parallel 1D FFT) as it offers the best speed-up with a higher number of threads. The worst algorithm seems to be the one which parallelizes only the 1D FFT (Serial 2D FFT w/ Parallel 1D FFT) since the speed-up seems to decrease when the number of threads is high.

## References
- [Steve Brunton - Image Compression and the FFT](https://www.youtube.com/watch?v=gGEBUdM0PVc)
- [Steve Brunton - Image Compression and the FFT (Examples in Python)](https://www.youtube.com/watch?v=uB3v6n8t2dQ)
