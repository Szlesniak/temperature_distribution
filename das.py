import matplotlib.pyplot as plt
import numpy as np
# https://pbc.biaman.pl/Content/362/PDF/wlasciwosci_termofizyczne_zywnosci_cz1.pdf
grid_size = 11
step = 1.5

temperatures = np.zeros((grid_size, grid_size))
center = grid_size // 2  # Środek siatki

for i in range(grid_size):
    for j in range(grid_size):
        r = np.sqrt(((i - center) * step) ** 2 + ((j - center) * step) ** 2)

        if r <= 1.5:  # Środek
            temperatures[i, j] = 26
        elif i == 0 or i == grid_size-1 or j == 0 or j == grid_size-1:
            temperatures[i, j] = 35
        elif r <= 4.5:  # Pierścień 3
            temperatures[i, j] = 28
        elif r <= 8:  # Pierścień 2
            temperatures[i, j] = 30
        else:
            temperatures[i, j] = 35


laplacian = np.zeros_like(temperatures)

for i in range(1, grid_size - 1):
    for j in range(1, grid_size - 1):
        laplacian[i, j] = (
            (temperatures[i + 1, j] - 2 * temperatures[i, j] + temperatures[i - 1, j]) / step ** 2 +
            (temperatures[i, j + 1] - 2 * temperatures[i, j] + temperatures[i, j - 1]) / step ** 2
        )


k = 1.15
q = -k * laplacian

print("Temperatury (°C):")
print(temperatures)
print("\nLaplasjan (°C/cm²):")

# Kontrolowanie liczby cyfr po przecinku
np.set_printoptions(precision=2, suppress=True)

print("Rozkład q(x, y) [W/m³]:")
print(q)
Q = q*10**4

stepm = 1.5 / 100

n = grid_size * grid_size
A = np.zeros((n,n))
b = np.zeros(n)

def idx(l, m):
    return l * grid_size + m

for i in range(grid_size):
    for j in range(grid_size):
        current = idx(i, j)

        if i == 0 or i == grid_size - 1 or j == 0 or j == grid_size - 1:
            A[current, current] = 1
            b[current] = 35
        else:
            A[current, current] = -4
            A[current, idx(i + 1, j)] = 1
            A[current, idx(i - 1, j)] = 1
            A[current, idx(i, j + 1)] = 1
            A[current, idx(i, j - 1)] = 1
            b[current] = -Q[i, j] * stepm ** 2 / k  # Uwzględnienie źródła ciepła q(x, y)
# print("A:\n")
# np.set_printoptions(precision=2, linewidth=140)
# print(A[20:50,20:50])
# print(b)

T = np.linalg.solve(A, b)
temperature_distribution = T.reshape((grid_size, grid_size))
print(temperatures)
print("Rozkład temperatury (°C):")
print(temperature_distribution)
x = np.linspace(0, 16, grid_size)
y = np.linspace(0, 16, grid_size)
X, Y = np.meshgrid(x, y)

plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, temperature_distribution, cmap='hot', levels=20)
plt.colorbar(contour, label='Temperatura [°C]')
plt.title('Rozkład temperatury w bigosie z uwzględnieniem q(x, y)')
plt.xlabel('X [cm]')
plt.ylabel('Y [cm]')
plt.show()

A1 = np.zeros((n,n))
b1 = np.zeros(n)
k1 = 1.15

for i in range(grid_size):
    for j in range(grid_size):
        current = idx(i, j)

        if i == 0 or i == grid_size - 1 or j == 0 or j == grid_size - 1:
            A1[current, current] = 1
            b1[current] = 49
        else:
            # Węzły wewnętrzne
            A1[current, current] = -4
            A1[current, idx(i + 1, j)] = 1
            A1[current, idx(i - 1, j)] = 1
            A1[current, idx(i, j + 1)] = 1
            A1[current, idx(i, j - 1)] = 1
            b1[current] = -Q[i, j] * 1.5 * stepm ** 2 / k1

T1 = np.linalg.solve(A1, b1)
temperature_distribution1 = T1.reshape((grid_size, grid_size))
print("Rozkład temperatury (°C):")
print(temperature_distribution1)
x = np.linspace(0, 16, grid_size)
y = np.linspace(0, 16, grid_size)
X, Y = np.meshgrid(x, y)

plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, temperature_distribution1, cmap='hot', levels=20)
plt.colorbar(contour, label='Temperatura [°C]')
plt.title('Rozkład temperatury w bigosie z uwzględnieniem q(x, y)')
plt.xlabel('X [cm]')
plt.ylabel('Y [cm]')
plt.show()

plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, temperatures, cmap='hot', levels=20)
plt.colorbar(contour, label='Temperatura [°C]')
plt.title('Zmierzone temperatury:')
plt.xlabel('X [cm]')
plt.ylabel('Y [cm]')
plt.show()