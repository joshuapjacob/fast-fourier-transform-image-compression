#include <complex>
#include <valarray>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <thread>
#include <execution>
#include <algorithm>
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;

// 1D COOLEY-TUKEY RADIX-2 DECIMATION-IN-FREQUENCY FAST FOURIER TRANSFORM ------

void fft_radix2(std::valarray<Complex> &x) {
    const size_t N = x.size();
    if (N <= 1) return;
    std::valarray<Complex> even = x[std::slice(0, N/2, 2)];
    std::valarray<Complex> odd = x[std::slice(1, N/2, 2)];
    fft_radix2(even);
    fft_radix2(odd);
    for (size_t k = 0; k < N/2; k++) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

void inverse_fft_radix2(std::valarray<Complex> &x) {
    x = x.apply(std::conj);
    fft_radix2(x);
    x = x.apply(std::conj);
    x /= x.size();
}

// OPTIMIZED PARALLEL 1D FAST FOURIER TRANSFORM --------------------------------

void parallel_fft_radix2(std::valarray<Complex> &x) {
    // TODO
    const size_t N = x.size();
    if (N <= 1) return;
    std::valarray<Complex> even = x[std::slice(0, N/2, 2)];
    std::valarray<Complex> odd = x[std::slice(1, N/2, 2)];
    fft_radix2(even);
    fft_radix2(odd);
    for (size_t k = 0; k < N/2; k++) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

void parallel_inverse_fft_radix2(std::valarray<Complex> &x) {
    x = x.apply(std::conj);
    parallel_fft_radix2(x);
    x = x.apply(std::conj);
    x /= x.size();
}

// 2D FAST FOURIER TRANSFORM UTILITY FUNCTIONS ---------------------------------

void fft_rows(
        std::valarray<Complex> &x, int width, int height, int start, int end
    ) {
    for (int i = start; i < end; i++) {
        std::valarray<Complex> sliced = x[std::slice(i*width,width,1)];
        fft_radix2(sliced);
        for (size_t k = 0; k < sliced.size(); k++) {
            x[i*width+k] = sliced[k];
        }
    }
}

void fft_columns(
        std::valarray<Complex> &x, int width, int height, int start, int end
    ) {
    for (int j = start; j < end; j++) {
        std::valarray<Complex> sliced = x[std::slice(j,height,width)];
        fft_radix2(sliced);
        for (size_t k = 0; k < sliced.size(); k++)
        {
            x[j+k*width] = sliced[k];
        }
    }
}

void inverse_fft_rows(
        std::valarray<Complex> &x, int width, int height, int start, int end
    ) {
    for (int i = start; i < end; i++) {
        std::valarray<Complex> sliced = x[std::slice(i*width,width,1)];
        inverse_fft_radix2(sliced);
        for (size_t k = 0; k < sliced.size(); k++) {
            x[i*width+k] = sliced[k];
        }
    }
}

void inverse_fft_columns(
        std::valarray<Complex> &x, int width, int height, int start, int end
    ) {
    for (int j = start; j < end; j++) {
        std::valarray<Complex> sliced = x[std::slice(j,height,width)];
        inverse_fft_radix2(sliced);
        for (size_t k = 0; k < sliced.size(); k++)
        {
            x[j+k*width] = sliced[k];
        }
    }
}

// SERIAL 2D FAST FOURIER TRANSFORM --------------------------------------------

void fft_image(std::valarray<Complex> &x, int width, int height) {
    fft_rows(x, width, height, 0, height);
    fft_columns(x, width, height, 0, width);
}

void inverse_fft_image(std::valarray<Complex> &x, int width, int height) {
    inverse_fft_columns(x, width, height, 0, width);
    inverse_fft_rows(x, width, height, 0, height);
}

// SIMPLE PARALLEL 2D FAST FOURIER TRANSFORM -----------------------------------

void simple_parallel_fft_image(
        std::valarray<Complex> &x, int width, int height, uint num_threads
    ) {
    uint block_size = height/num_threads;
    std::vector<std::thread> threads(num_threads - 1);
    uint block_start = 0;
    uint block_end = block_size;
    for (uint i = 0; i < (num_threads - 1); i++) {
        threads[i] = std::thread(
            fft_rows,
            std::ref(x), width, height,
            block_start, block_end
        );
        block_start += block_size;
        block_end += block_size;
    }
    fft_rows(x, width, height, block_start, height);
    std::for_each(threads.begin(), threads.end(),
        std::mem_fn(&std::thread::join)
    );
    block_size = width/num_threads;
    block_start = 0;
    block_end = block_size;
    for (uint i = 0; i < (num_threads - 1); i++) {
        threads[i] = std::thread(
            fft_columns,
            std::ref(x), width, height,
            block_start, block_end
        );
        block_start += block_size;
        block_end += block_size;
    }
    fft_columns(x, width, height, block_start, width);
    std::for_each(threads.begin(), threads.end(),
        std::mem_fn(&std::thread::join)
    );
}

void simple_parallel_inverse_fft_image(
        std::valarray<Complex> &x, int width, int height, uint num_threads
    ) {
    uint block_size = width/num_threads;
    std::vector<std::thread> threads(num_threads - 1);
    uint block_start = 0;
    uint block_end = block_size;
    for (uint i = 0; i < (num_threads - 1); i++) {
        threads[i] = std::thread(
            inverse_fft_columns,
            std::ref(x), width, height,
            block_start, block_end
        );
        block_start += block_size;
        block_end += block_size;
    }
    inverse_fft_columns(x, width, height, block_start, width);
    std::for_each(threads.begin(), threads.end(),
        std::mem_fn(&std::thread::join)
    );
    block_size = height/num_threads;
    block_start = 0;
    block_end = block_size;
    for (uint i = 0; i < (num_threads - 1); i++) {
        threads[i] = std::thread(
            inverse_fft_rows,
            std::ref(x), width, height,
            block_start, block_end
        );
        block_start += block_size;
        block_end += block_size;
    }
    inverse_fft_rows(x, width, height, block_start, height);
    std::for_each(threads.begin(), threads.end(),
        std::mem_fn(&std::thread::join)
    );

}

// OPTIMIZED PARALLEL 2D FAST FOURIER TRANSFORM --------------------------------

// TODO

// COMPRESSION -----------------------------------------------------------------

void quick_compress_block(
        std::valarray<Complex> &x, double cut_off, uint start, uint end
    ) {
    for (size_t k = start; k < end; k++) {
        if (abs(x[k]) < cut_off) x[k] = 0;
    }
}

void quick_compress(
        std::valarray<Complex> &x, double cut_off, uint num_threads
    ) {
    uint block_size = x.size()/num_threads;
    std::vector<std::thread> threads(num_threads - 1);
    uint block_start = 0;
    uint block_end = block_size;
    for (uint i = 0; i < (num_threads - 1); i++) {
        threads[i] = std::thread(
            quick_compress_block,
            std::ref(x), cut_off,
            block_start, block_end
        );
        block_start += block_size;
        block_end += block_size;
    }
    quick_compress_block(x, cut_off, block_start, x.size());
    std::for_each(threads.begin(), threads.end(),
        std::mem_fn(&std::thread::join)
    );
}

double find_cut_off(std::valarray<Complex> &x, double keep) {
    // I did not try to optimize this function.
    std::valarray<Complex> *sorted_x = new std::valarray(x);
    struct {
        bool operator()(Complex a, Complex b) const {
            return abs(a) > abs(b);
        }
    } comp;
    std::sort(
        std::execution::par_unseq,
        std::begin(*sorted_x), std::end(*sorted_x),
        comp
    );
    size_t cut_off_index = std::round(keep*x.size());
    return abs((*sorted_x)[cut_off_index]);
}



// UTILITY FUNCTIONS -----------------------------------------------------------

template<typename F, typename... Args>
double time_it(F func, Args&&... args){
    auto start_time = std::chrono::high_resolution_clock::now();
    func(std::forward<Args>(args)...);
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time
    );
    return duration.count() / 1000.;
}

std::valarray<Complex>* read_img(
        const char *filename, int &width, int &height
    ) {
    int channels;
    uint8_t* rgb_img = stbi_load(filename, &width, &height, &channels, 3);

    std::cout << "Image Size: " << width << "x" << height <<
        " (" << channels << " channels)" << std::endl;

    size_t rgb_img_size = width*height*channels;
    size_t grey_img_size = width*height;

    Complex *img_array =  new Complex[grey_img_size];

    size_t k = 0;
    // TODO: Parallelize.
    for (uint8_t* p = rgb_img; p != rgb_img + rgb_img_size; p += channels) {
        img_array[k] = 0.2989*(*p) + 0.5870*(*(p + 1)) + 0.1140*(*(p + 2));
        k++;
    }

    std::valarray<Complex>* img_data = new std::valarray<Complex>(
        img_array, grey_img_size
    );

    delete [] img_array;
    stbi_image_free(rgb_img);

    return img_data;
}

void write_img(
        std::valarray<Complex> &img_data, int width, int height,
        std::string filename
    ) {
    size_t out_img_size = width*height;
    uint8_t out_img[out_img_size];
    // TODO: Parallelize.
    for (size_t k = 0; k < out_img_size; k++) {
        out_img[k] = (uint8_t) std::max(std::min(255., img_data[k].real()),0.);
    }
    const char *c = filename.c_str();
    stbi_write_png(c, width, height, 1, out_img, width);

}

// MAIN ------------------------------------------------------------------------

int main(int argc, char const *argv[]) {

    if (argc != 3) {
        std::cerr <<
            "Please provide one input image filename and "
            "fraction of data to keep (between 0 and 1)."
        << std::endl;
        exit(1);
    }

    const char *input_image_filename = argv[1];
    if (not std::filesystem::exists(input_image_filename)) {
        std::cerr << "Image file not found." << std::endl;
        exit(1);
    }

    const double keep = atof(argv[2]);
    if (keep > 1.0 || keep < 0.0) {
        std::cerr <<
            "Fraction of data to keep should be between 0 and 1."
        << std::endl;
        exit(1);
    }

    uint const n = 4; // std::thread::hardware_concurrency();
    std::cout << "Number of threads: " << n << std::endl << std::endl;

    // TEST: Serial 2D FFT w/ Serial 1D FFT

    int width, height;
    std::valarray<Complex> *img_data = read_img(
        input_image_filename, width, height
    );

    std::cout << "FFT Image Duration: " <<
        time_it(fft_image, *img_data, width, height) << "s" << std::endl;

    double cut_off = find_cut_off(*img_data, keep); // WARNING: SLOW

    std::cout << "Quick Compression Duration: " <<
        time_it(quick_compress, *img_data, cut_off, n) << "s" << std::endl;

    std::cout << "Inverse FFT Image Duration: " <<
        time_it(inverse_fft_image, *img_data, width, height) << "s" << std::endl;

    write_img(*img_data, width, height, "output1.png");
    delete img_data;

    std::cout << std::endl;

    // TEST: Simple Parallel 2D FFT w/ Serial 1D FFT

    img_data = read_img(
        input_image_filename, width, height
    );

    std::cout << "Naive Parallel FFT Image Duration: " <<
        time_it(simple_parallel_fft_image, *img_data, width, height, n) << "s" << std::endl;

    std::cout << "Quick Compression Duration: " <<
        time_it(quick_compress, *img_data, cut_off, n) << "s" << std::endl;

    std::cout << "Naive Parallel Inverse FFT Image Duration: " <<
        time_it(simple_parallel_inverse_fft_image, *img_data, width, height, n) << "s" << std::endl;

    write_img(*img_data, width, height, "output2.png");

    // TEST: Optimized Parallel 2D FFT w/ Serial 1D FFT
    // TODO

    // TEST: Serial 2D FFT w/ Parallel 1D FFT
    // TODO

    // TEST: Simple Parallel 2D FFT w/ Parallel 1D FFT
    // TODO

    // TEST: Optimized Parallel 2D FFT w/ Parallel 1D FFT
    // TODO

    delete img_data;
    return 0;
}
