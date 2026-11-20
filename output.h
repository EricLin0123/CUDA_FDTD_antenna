#ifndef FDTD_OUTPUT_H
#define FDTD_OUTPUT_H
#include <cstdio>
#include <cstdlib>
#include <string>

inline FILE *checked_open(const char *path, const char *mode) {
    FILE *file = fopen(path, mode);
    if (!file) { perror(path); exit(EXIT_FAILURE); }
    return file;
}

// NumPy v1.0: native-endian float64, C order, (frame, x, y).
class FieldOutput {
    FILE *file;
public:
    FieldOutput(const char *path, int frames, int nx, int ny) {
        file = checked_open(path, "wb");
        unsigned short endian = 1;
        std::string header = "{'descr': '";
        header += (*reinterpret_cast<unsigned char *>(&endian) ? "<f8" : ">f8");
        header += "', 'fortran_order': False, 'shape': (" + std::to_string(frames)
            + ", " + std::to_string(nx) + ", " + std::to_string(ny) + "), }";
        header.append((64 - ((10 + header.size() + 1) % 64)) % 64, ' ');
        header += '\n';
        unsigned char prefix[] = {0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0,
            static_cast<unsigned char>(header.size() & 255),
            static_cast<unsigned char>(header.size() >> 8)};
        if (fwrite(prefix, 1, sizeof(prefix), file) != sizeof(prefix)
            || fwrite(header.data(), 1, header.size(), file) != header.size()) {
            perror("writing NumPy header"); exit(EXIT_FAILURE);
        }
    }
    void write(const double *values, size_t count) {
        if (fwrite(values, sizeof(double), count, file) != count) {
            perror("writing fields"); exit(EXIT_FAILURE);
        }
    }
    void close() {
        if (fclose(file)) { perror("closing fields"); exit(EXIT_FAILURE); }
    }
};
#endif
