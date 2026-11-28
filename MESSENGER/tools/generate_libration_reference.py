"""Независимые эталоны G20q: интеграл по M, Ньютон для E, квадратура Гаусса.
Рабочая C++ модель использует интеграл по E и составную формулу Симпсона.
Для пересчёта эталонов требуется NumPy; при сборке C++ Python не нужен.
"""
from pathlib import Path
import csv
import numpy as np

root = Path(__file__).resolve().parents[1]
nodes, weights = np.polynomial.legendre.leggauss(128)
M = np.pi * (nodes + 1)
with (root / 'tests/libration_reference.csv').open('w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)
    writer.writerow(['eccentricity', 'q', 'G20q'])
    for e in (0., .1, .2056317, .3):
        E = M.copy()
        for _ in range(20):
            correction = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
            E -= correction
            if np.max(np.abs(correction)) < 1e-15:
                break
        f = np.arctan2(np.sqrt(1-e*e)*np.sin(E), np.cos(E)-e)
        for q in range(-4, 7):
            value = .5 * np.sum(weights * np.cos(2*f-(2+q)*M) / (1-e*np.cos(E))**3)
            writer.writerow([format(e, '.17g'), q, format(value, '.17g')])
print('44 независимых эталона записаны')
