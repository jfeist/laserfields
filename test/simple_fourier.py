import sys
import numpy as np

ts, Et, *_ = np.loadtxt(sys.stdin.readlines(), unpack=True)

# we want this frequency resolution
dw_target = 0.0015197108631458
dt = ts[1] - ts[0]
nn = int(2*np.pi/(dw_target*dt))
ws = np.fft.fftfreq(nn, dt) * 2 * np.pi
Ef = np.fft.fft(Et, nn)
# this scaling ensures that the functions are normalized correctly, and that the 
# phase is correct for the time interval starting at ts[0]
Ef *= np.exp(-1j*ws*ts[0]) * np.sqrt(dt**2 / (2*np.pi))
# assert np.allclose(np.trapezoid(Et**2,ts), np.trapezoid(abs(Ef)**2, dx=ws[1]-ws[0]))

inds = (ws>0) & (ws < 300/27.211386245988)
ws, Ef = ws[inds], Ef[inds]

np.savetxt(sys.stdout.buffer, np.column_stack((ws, Ef.real, Ef.imag)),
           header="omega                  |E(w)|^2           arg(E(w))",
           fmt="%.14e %.14e %+.14e")
