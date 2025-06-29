#include "utils.h"

// Función para generar cadenas aleatorias
std::string generateRandomString(const std::string& alphabet, int length) {
    std::string result;
    result.reserve(length);

    std::mt19937 rng(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_int_distribution<int> dist(0, alphabet.length() - 1);

    for (int i = 0; i < length; ++i) {
        result.push_back(alphabet[dist(rng)]);
    }
    return result;
} 