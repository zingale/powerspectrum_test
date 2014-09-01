# create some "fake" 1-d multimode data and create a power spectrum to
# illustrate that we have the scaling and normalization correct

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

# number of discrete points in real space to sample
N = 64

index = -4.0

xmin = 0.0
xmax = 10.0

L = xmax - xmin

phi = np.zeros((N), dtype=np.float64)

x = (np.arange(N)+0.5)*(xmax - xmin)/N + xmin

modes = 16

A_0 = 1000.0


# amplitude of smallest k
print "min{k_x} = ", 1/L

A_ksmall = A_0*(1/L)**index
print "amplitude of smallest k = ", np.sqrt(A_ksmall)
print "power at smallest k =     ", A_ksmall


# compute the function we will find the power spectrum of
for m in range(1,modes+1):
    k_x = m/L

    A = np.sqrt(A_0*k_x**index)
            
    phi += A*np.sin(2.0*np.pi*k_x*x)


#---------------------------------------------------------------------------
# plot the real-space function
#---------------------------------------------------------------------------

plt.plot(x, phi)
plt.xlabel(r"$x$")
plt.ylabel(r"$f(x)$")

plt.tight_layout()
plt.savefig("phi_1d.png")
plt.savefig("phi_1d.eps")


#---------------------------------------------------------------------------
# now do the Fourier transform
#---------------------------------------------------------------------------
phi_hat = np.fft.fft(phi)


# Parseval's theorem: sum of |phi(x)|**2 = sum of (1/N) |phi_hat(k)|**2
parseval_thm_realspace = np.sum(np.abs(phi)**2)
parseval_thm_fourier = np.sum(np.abs(phi_hat)**2)/N

print "sum of |phi(x)|**2      = {}".format(parseval_thm_realspace)
print "sum of |phihat(k)|**2/N = {}".format(parseval_thm_fourier)
print " "
print "ratio = {}".format(parseval_thm_realspace/parseval_thm_fourier)
print "DC Offset = {}".format(phi_hat[0])


#---------------------------------------------------------------------------
# normalization for power spectrum later
#---------------------------------------------------------------------------

# only positive is really unique
phi_hat = 2.0*phi_hat[0:N/2+1]

phi_hat = phi_hat/N


#---------------------------------------------------------------------------
# get the wavenumbers in physical [cm^{-1}] units
#---------------------------------------------------------------------------
kx = np.fft.rfftfreq(N)
kx = kx*N/L


#---------------------------------------------------------------------------
# plot the transformed function
#---------------------------------------------------------------------------

plt.clf()

plt.plot(kx, np.abs(phi_hat))
plt.xlabel(r"$k$")
plt.ylabel(r"$|\mathcal{F}(k)|$")

plt.tight_layout()
plt.savefig("phihat_1d.png")
plt.savefig("phihat_1d.eps")


#---------------------------------------------------------------------------
# plot the power spectrum
#---------------------------------------------------------------------------

plt.clf()

E_spectrum = np.abs(phi_hat)**2
plt.loglog(kx, E_spectrum)

ii = np.argmax(E_spectrum)
kmax = kx[ii]
Emax = E_spectrum[ii]

print "maximum E = {} at k = {}".format(Emax, kmax)

plt.loglog(kx, Emax*(kx/kmax)**index, ls=":", color="0.5")
plt.ylim(1.e-10*Emax, 1.1*Emax)

plt.xlabel("k")
plt.ylabel("E(k)dk")

plt.tight_layout()
plt.savefig("ps1d.png")
plt.savefig("ps1d.eps")


