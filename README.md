# fast-fourier-transform
CSE305 Final Project

N.B. The program does not actually output a compressed jpeg image.

## Compilation & Execution
```
$ g++ -O3 main.cpp -Wall -Wextra -pedantic -lpthread -D_GLIBCXX_PARALLEL -fopenmp -lttb -std=c++17 -o main.out
$ ./main.out images/philip.png 0.05 # Keeps only 5% of the data in philip.png
```
Output image is stored as `output.png`.

Image width and height must be a power of two.

"usage: ./a.out <size of array> <number of threads>"
## TODO
- Do not slice in 2d transform
- Perform  a  detailed  comparison  of  the  resulting  algorithms.
- Plot time vs. threads
- Plot speed-up vs. threads
- show different compression rates for different images

## Bonus: Code Quality
- No compiler warnings (with `-Wall -Wextra -pedantic`)
- No memory leaks