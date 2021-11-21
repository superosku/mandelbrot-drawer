#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include "lodepng.h"
#include <stdlib.h>

class Complex {
public:
    double a;
    double b;

    Complex() {
        a = 0.0;
        b = 0.0;
    }

    Complex(double a_, double b_) {
        a = a_;
        b = b_;
    }

    Complex operator*(const Complex& other) {
        return Complex(
            a * other.a - b * other.b,
            a * b + b * a
        );
    }

    Complex operator+(const Complex& other) {
        return Complex(
            a + other.a,
            b + other.b
        );
    }

    double abs() {
        return sqrt(a * a + b * b);
    }
};


std::ostream& operator<<(std::ostream& os, const Complex& c)
{
    os << "(" << c.a << ", " << c.b << ")";
    return os;
}


class MandelbrotImage {
    int res;
    std::vector<double> data;

    Complex center;

public:
    double zoom;

    MandelbrotImage(int res_, Complex center_, double zoom_) {
        res = res_;
        center = center_;
        zoom = zoom_;
        data.resize(res * res);
        for (int i = 0; i < res * res; i ++) {
            data[i] = 0;
        }
    }

    void calculate() {
        for (int i = 0; i < res * res; i ++) {
            data[i] = 0;
        }

        std::vector<Complex> c_data;
        std::vector<Complex> z_data;
        c_data.resize(res * res);
        z_data.resize(res * res);

        for (int j = 0; j < res; j ++) {
            for (int i = 0; i < res; i ++) {
                c_data[i + j * res] = Complex(
                    zoom * 2.0 * i / (res - 1) - zoom + center.a,
                    zoom * 2.0 * j / (res - 1) - zoom + center.b
                );
                z_data[i + j * res] = Complex();
            }
        }

        int total_change_count = 0;
        int k = 0;
        for (; k < 2000; k ++) {
            int loop_change_count = 0;
            for (int j = 0; j < res; j ++) {
                for (int i = 0; i < res; i ++) {
                    if (data[i + j * res] > 0) {
                        continue;
                    }

                    z_data[i + j * res] = z_data[i + j * res] * z_data[i + j * res] + c_data[i + j * res];

                    if (z_data[i + j * res].abs() > 2) {
                        data[i + j * res] = k;
                        total_change_count += 1;
                        loop_change_count += 1;
                    }
                }
            }
            // Stop if we have written most of the image and are not writing any more
            if (
                total_change_count > res * res / 25 &&
                (
                    loop_change_count < res * res / 10000 ||
                    loop_change_count < 10
                )
            ) {
//                std::cout << "thing: " << loop_change_count << " " << k << std::endl;
                break;
            }
        }
        for (unsigned long i = 0; i < data.size(); i ++) {
            if (data[i] == 0) {
                data[i] = k + 1;
            }
        }
        double max_val = 0;
        double min_val = 99999999;
        for (unsigned long i = 0; i < data.size(); i ++) {
            max_val = data[i] > max_val ? data[i] : max_val;
            min_val = data[i] < min_val ? data[i] : min_val;
        }
        for (unsigned long i = 0; i < data.size(); i ++) {
            data[i] = (data[i] - min_val) / (max_val - min_val);
        }
    }

    void save(int batch, int id, int aliasing) {
        // Output image resolution is res / aliasing
        std::vector<unsigned char> image;
        int output_res = res / aliasing;
        image.resize(output_res * output_res * 4);

        for (int i = 0; i < output_res; i ++) {
            for (int j = 0; j < output_res; j ++) {
                double val = 0.0;
                // Loop anti aliasing filter
                for (int k = 0; k < aliasing; k ++) {
                    for (int l = 0; l < aliasing; l ++) {
                        val += data[
                            (i * aliasing + k) +
                            (j * aliasing + l) * res
                        ];
                    }
                }
                val = val / (aliasing * aliasing);
                unsigned char cval = 255.0 * val;
                image[(i + j * output_res) * 4 + 0] = cval;
                image[(i + j * output_res) * 4 + 1] = cval;
                image[(i + j * output_res) * 4 + 2] = cval;
                image[(i + j * output_res) * 4 + 3] = 255;
            }
        }

//        for (unsigned long i = 0; i < data.size(); i ++) {
//            unsigned char val = 255.0 * data[i];
//            image[i * 4 + 0] = val;
//            image[i * 4 + 1] = val;
//            image[i * 4 + 2] = val;
//            image[i * 4 + 3] = 255;
//        }

        unsigned error = lodepng::encode(
            (
                "cppoutput/image-" +
                std::to_string(batch) +
                "-" +
                std::to_string(id) +
                ".png"
            ),
            image,
            output_res,
            output_res
        );

        if (error) {
            std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
        }
    }

    MandelbrotImage get_zoomed_img() {
        MandelbrotImage best(32, center, zoom * 0.75);
        best.calculate();
        int best_score = best.unique_count();

        for (int i = 0; i < 10; i ++) {
            double randa = (((double)rand()) / RAND_MAX - 0.5) * zoom * 0.75;
            double randb = (((double)rand()) / RAND_MAX - 0.5) * zoom * 0.75;
//            std::cout << randa << " " << randb << std::endl;
            MandelbrotImage candidate(
                32,
                Complex(
                    center.a + randa,
                    center.b + randb
                ),
                zoom * 0.75
            );
            candidate.calculate();
            int candidate_score = candidate.unique_count();
            if (candidate_score > best_score) {
                best_score = candidate_score;
                best = candidate;
            }
        }

        return MandelbrotImage(
            res,
            best.center,
            best.zoom
        );
    }

    int unique_count() {
        std::map<double, bool> m;

        for (unsigned long i = 0; i < data.size(); i ++) {
            m[data[i]] = true;
        }

        return m.size();
    }
};


int main() {
    std::cout << "Hello world!" << std::endl;



    #pragma omp parallel for schedule(dynamic) num_threads(8)
    for (int b = 0; b < 16; b ++) {
        MandelbrotImage img(2 * 4 * 256, Complex(-0.75, 0.0), 1.5);
        for (int i = 0; i < 100; i ++) {
//            img.calculate();
//            img.save(b, i, 4);
            img = img.get_zoomed_img();
        }
        std::cout << "got the coords" << std::endl;

        img.zoom = 1.5;
        for (int i = 0; i < 100; i ++) {
            img.calculate();
            img.save(b, i, 4);
            img.zoom = img.zoom * 0.75;
            std::cout << "img written" << std::endl;
        }
    }

    return 0;
}

