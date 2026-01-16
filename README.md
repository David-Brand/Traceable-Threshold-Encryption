# Traceable-Threshold-Encryption

An open-source implementation of a **Traceable Threshold Key Encapsulation Mechanism (TTT-KEM)**, which were introduced by Boneh, Partap, and Rotem [CRYPTO'24].

⚠️ **Disclaimer**:  
This code is provided for research and educational purposes only. It has not been formally verified or security-audited and should not be used in production systems.

---

## Features

- Traceable threshold key encapsulation
- Pairing-based cryptography backend
- Modular structure for easy modification and extension
- Suitable as a reference or prototype implementation

---

## Dependencies

This project depends on the following external library:

### mcl

- **mcl** (Multiprecision Cryptography Library)
- Used for elliptic curve and pairing-based cryptographic operations

Repository: https://github.com/herumi/mcl

Please make sure `mcl` is available before building this project.

---

## Building the Project

### Prerequisites

- C++ compiler with **C++17** support
- `cmake`
- `mcl` library

---

### build

- To build the project, set the path to the `mcl` base directory in the Cmakelists.txt file.
- The `mcl` library needs to be built specifically for the compiler and system you are using (I used msvc 19), if this is not met, the project won't build.

Then, create a `build` directory, navigate into it, and build the project using cmake.

```bash
mkdir build && cd build
cmake ..
cmake --build . --config Release
```

---

## License

This project is licensed under the MIT License.
See the LICENSE.txt file for details.
