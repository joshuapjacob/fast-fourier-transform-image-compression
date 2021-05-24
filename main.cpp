#include <complex>
#include <valarray>
#include <iostream>
#include <chrono>
#include <filesystem>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;

// COOLEY-TUKEY FAST FOURIER TRANSFORM RADIX-2 ---------------------------------

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

// TODO: IMAGE DATA MATRIX -----------------------------------------------------
// (https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op)

// template<typename ty>
// class Matrix {
//     std::valarray<ty> data;
//     int dim;
//  public:
//     Matrix(int r, int c) : data(r*c), dim(c) {}
//     ty& operator()(int r, int c) {return data[r*dim + c];}
// };

// 2D FAST FOURIER TRANSFORM ---------------------------------------------------

void fft_image(std::valarray<Complex> &x, int width, int height) {
    // TODO: FIX!

    for (size_t i = 0; i < height; i++) {
        std::valarray<Complex> sliced = x[std::slice(i*width,width,1)];
        fft_radix2(sliced);
        for (size_t k = 0; k < sliced.size(); k++) {
            x[i*width+k] = sliced[k];
        }
    }

    for (size_t j = 0; j < width; j++) {
        std::valarray<Complex> sliced = x[std::slice(j,height,width)];
        fft_radix2(sliced);
        for (size_t k = 0; k < sliced.size(); k++)
        {
            x[j+k*width] = sliced[k];
        }
    }

    for (size_t i = 0; i < 10; i++)
    {
        std::cout << x[i] << std::endl;
    }
}

void inverse_fft_image(std::valarray<Complex> &x, int width, int height) {
    // TODO: FIX!

    // for (size_t j = 0; j < width; j++) {
    //     std::valarray<Complex> sliced = x[std::slice(j,height,width)];
    //     inverse_fft_radix2(sliced);
    //     for (size_t k = 0; k < sliced.size(); k++)
    //     {
    //         x[j+k*width] = sliced[k];
    //     }
    // }

    // for (size_t i = 0; i < height/2; i++) {
    //     std::valarray<Complex> sliced = x[std::slice(i*width,width,1)];
    //     inverse_fft_radix2(sliced);
    //     for (size_t k = 0; k < sliced.size(); k++) {
    //         x[i*width+k] = sliced[k];
    //     }
    // }
}

// COMPRESSION -----------------------------------------------------------------

void compress(std::valarray<Complex> &x) {
    // TODO: FIX!
    for (size_t k = 0; k < x.size(); k++) {
        // std::cout << x[k] << " ABS: " << abs(x[k]) << std::endl;
        if (abs(x[k]) < 100000) {
            x[k] = 0;
        }
    }
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

std::valarray<Complex>* read_img(const char *filename, int &width, int &height) {
    int channels;
    uint8_t* rgb_img = stbi_load(filename, &width, &height, &channels, 3);

    std::cout << "Image Size: " << width << "x" << height <<
        " (" << channels << " channels)" << std::endl;

    size_t rgb_img_size = width*height*channels;
    size_t grey_img_size = width*height;

    Complex *img_array =  new Complex[grey_img_size];

    size_t k = 0;
    for (uint8_t* p = rgb_img; p != rgb_img + rgb_img_size; p += channels) {
        img_array[k] = (*p + *(p + 1) + *(p + 2))/3.0; //TODO: Use RGB Weights! rgb_weights = [0.2989, 0.5870, 0.1140]?
        k++;
    }

    std::valarray<Complex>* img_data = new std::valarray<Complex>(
        img_array, grey_img_size
    );

    delete [] img_array;
    stbi_image_free(rgb_img);

    return img_data;
}

void write_img(std::valarray<Complex> &img_data, int width, int height) {
    size_t out_img_size = width*height;
    uint8_t out_img[out_img_size];
    for (size_t k = 0; k < out_img_size; k++) {
        out_img[k] = (uint8_t)std::round(img_data[k].real());
    }
    stbi_write_png("output.png", width, height, 1, out_img, width);

}

// MAIN ------------------------------------------------------------------------

int main(int argc, char const *argv[]) {

    if (argc =! 2) {
        std::cout << "Please provide one input image filename." << std::endl;
        exit(1);
    }

    const char *input_image_filename = argv[1];
    if (not std::filesystem::exists(input_image_filename)) {
        std::cout << "Image file not found." << std::endl;
        exit(1);
    }

    int width, height;
    std::valarray<Complex> *img_data = read_img(
        input_image_filename, width, height
    );

    std::cout << "FFT Image Duration: " <<
        time_it(fft_image, *img_data, width, height) << "s" << std::endl;

    // std::cout << "Compression Duration:" <<
    //     time_it(compress, *img_data) << "s" << std::endl;

    std::cout << "Inverse FFT Image Duration: " <<
        time_it(inverse_fft_image, *img_data, width, height) << "s" << std::endl;

    // std::cout << "FFT Radix-2 Duration: " <<
    //     time_it(fft_radix2, *img_data) << "s" << std::endl;
    // std::cout << "Inverse FFT Radix-2 Duration: " <<
    //     time_it(inverse_fft_radix2, *img_data) << "s" << std::endl;

    write_img(*img_data, width, height);
    delete img_data;

    return 0;
}
