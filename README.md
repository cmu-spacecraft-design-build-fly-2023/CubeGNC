# CubeGNC

High-fidelity simulator for the development of GNC algorithms.

## Installation

Clone the repository and run in the root folder:
```bash
pip install .
```
In editable mode:

```bash
pip install -e .
```
## Documentation

See tests/test_spacecraft.py for a simple example

## TODO
- Sensor models: Magnetometer, Gyro, Camera outputs
- Magnetorquer
- Atmospheric drag - switch to nrlmsise00 for the density model
- Unit tests