# create some "fake" 3-d multimode data and create a power spectrum to
# illustrate that we have the scaling and normalization correct

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

# number of discrete points in real space to sample
N = 64

index = -4.0

xmin = ymin = zmin = 0.0
xmax = ymax = zmax = 10.0

L = np.array([xmax - xmin, ymax - ymin, zmax - zmin])

phi = np.zeros((N,N,N), dtype=np.float64)

x = (np.arange(N)+0.5)*(xmax - xmin)/N + xmin
y = (np.arange(N)+0.5)*(ymax - ymin)/N + ymin
z = (np.arange(N)+0.5)*(zmax - zmin)/N + zmin

x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

modes = 16

A_0 = 1000.0

# the minimum mode to consider -- to consider a single mode, set this
# to modes
mode_min = 1

# first pass -- count how many hits each k gets
weights = {}

for m in range(mode_min,modes+1):
    for n in range(mode_min,modes+1):
        for p in range(mode_min,modes+1):

            ii2 = m**2 + n**2 + p**2
            
            if not ii2 in weights.keys():
                weights[ii2] = 1
            else:
                weights[ii2] += 1


# amplitude of smallest k
print "min{k_x} = ", mode_min/L[0]

ii_min = 3*mode_min**2 
A_ksmall = A_0*np.sqrt((mode_min/L[0])**2 + 
                       (mode_min/L[1])**2 + 
                       (mode_min/L[2])**2)**index/weights[ii_min]

print "amplitude of smallest k = ", np.sqrt(A_ksmall)
print "power at smallest k =     ", A_ksmall
print "weights at smallest k =   ", weights[ii_min]

# compute the function we will find the power spectrum of
for m in range(mode_min,modes+1):
    k_x = m/L[0]
    print "m = ", m

    for n in range(mode_min,modes+1):
        k_y = n/L[1]

        for p in range(mode_min,modes+1):
            k_z = p/L[2]

            ii2 = m**2 + n**2 + p**2

            k = np.sqrt(k_x**2 + k_y**2 + k_z**2)
            A = np.sqrt(A_0*k**index/weights[ii2])
            print "({}, {}, {})".format(m,n,p), "k = ", k, "A = ", A
            phi += A*np.sin(2.0*np.pi*k_x*x3d + 
                            2.0*np.pi*k_y*y3d + 
                            2.0*np.pi*k_z*z3d)


#---------------------------------------------------------------------------
# plot the real-space function
#---------------------------------------------------------------------------
F = plt.figure()

grid = AxesGrid(F, (0.05, 0.1, 0.85, 0.85), 
                nrows_ncols = (1, 3),
                axes_pad = 1.1,
                label_mode = "all",
                share_all = False,
                cbar_location = "right",
                cbar_mode = "each",
                cbar_size = "3%",
                cbar_pad = "0%")


im = grid[0].imshow(phi[:,:,N/2], interpolation="nearest", origin="lower")
grid.cbar_axes[0].colorbar(im)

im = grid[1].imshow(phi[:,N/2,:], interpolation="nearest", origin="lower")
grid.cbar_axes[1].colorbar(im)

im = grid[2].imshow(phi[N/2,:,:], interpolation="nearest", origin="lower")
grid.cbar_axes[2].colorbar(im)

plt.savefig("phi.png")


#---------------------------------------------------------------------------
# now do the Fourier transform
#---------------------------------------------------------------------------
phi_hat = np.fft.fftn(phi)


# Parseval's theorem: sum of |phi(x)|**2 = sum of (1/N**3) |phi_hat(k)|**2
parseval_thm_realspace = np.sum(np.abs(phi)**2)
parseval_thm_fourier = np.sum(np.abs(phi_hat)**2)/N**3

print "sum of |phi(x)|**2         = {}".format(parseval_thm_realspace)
print "sum of |phihat(k)|**2/N**3 = {}".format(parseval_thm_fourier)
print " "
print "ratio = {}".format(parseval_thm_realspace/parseval_thm_fourier)
print "DC Offset = {}".format(phi_hat[0,0,0])


#---------------------------------------------------------------------------
# plot the Fourier transform
#---------------------------------------------------------------------------
plt.clf()

F = plt.figure()

grid = AxesGrid(F, (0.05, 0.1, 0.85, 0.85), 
                nrows_ncols = (1, 3),
                axes_pad = 1.1,
                label_mode = "all",
                share_all = False,
                cbar_location = "right",
                cbar_mode = "each",
                cbar_size = "3%",
                cbar_pad = "0%")


im = grid[0].imshow(np.abs(phi_hat[:,:,0]).T, interpolation="nearest", 
                    origin="lower")
grid[0].axes.xaxis.set_label_text("x")
grid[0].axes.yaxis.set_label_text("y")
grid.cbar_axes[0].colorbar(im)

im = grid[1].imshow(np.abs(phi_hat[:,0,:]).T, interpolation="nearest",
                    origin="lower")
grid[1].axes.xaxis.set_label_text("x")
grid[1].axes.yaxis.set_label_text("z")
grid.cbar_axes[1].colorbar(im)

im = grid[2].imshow(np.abs(phi_hat[0,:,:]).T, interpolation="nearest",
                    origin="lower")
grid[2].axes.xaxis.set_label_text("y")
grid[2].axes.yaxis.set_label_text("z")
grid.cbar_axes[2].colorbar(im)

plt.savefig("phihat.png")


#---------------------------------------------------------------------------
# normalization for power spectrum later
#---------------------------------------------------------------------------

# only a single octant is really unique
phi_hat = 8.0*phi_hat[0:N/2+1,0:N/2+1,0:N/2+1]

phi_hat = phi_hat/N**3


#---------------------------------------------------------------------------
# get the wavenumbers in physical [cm^{-1}] units
#---------------------------------------------------------------------------
kx = np.fft.rfftfreq(N)
kx = kx*N/L[0]

ky = np.fft.rfftfreq(N)
ky = ky*N/L[1]

kz = np.fft.rfftfreq(N)
kz = kz*N/L[2]

#---------------------------------------------------------------------------
# bin up |phi_hat| in terms of |k|
#---------------------------------------------------------------------------

dk = kx[1]-kx[0]

# we don't care about a wavenumber of 0, but we want the wavenumber
# bins centered on the physical values we care about
kmin = np.sqrt(kx[1]**2 + ky[1]**2 + kz[1]**2)-0.5*dk

num = int(np.sqrt(3)*N)
kmax = num*dk + kmin

bins = np.linspace(kmin, kmax, num+1, endpoint=True)

kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")

k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

whichbin = np.digitize(k.flat, bins)

ncount = np.bincount(whichbin)

E_spectrum = np.zeros(len(ncount)-1, dtype=np.float64)


# unresolved:  should I normalize by ncount[n]?  

for n in range(len(ncount)):
    if not ncount[n] == 0: 
        #E_spectrum[n-1] = np.sum((np.abs(phi_hat)).flat[whichbin==n]) /ncount[n]
        E_spectrum[n-1] = np.sum((np.abs(phi_hat)**2).flat[whichbin==n]) /ncount[n]


#---------------------------------------------------------------------------
# plot the power spectrum
#---------------------------------------------------------------------------

plt.clf()

k = 0.5*(bins[0:num] + bins[1:num+1])
k = k[0:len(E_spectrum)]

print k.shape
print E_spectrum.shape

print k[0], k[1], dk

plt.loglog(k, E_spectrum)

ii = np.argmax(E_spectrum)
kmax = k[ii]
Emax = E_spectrum[ii]

plt.loglog(k, Emax*(k/kmax)**index, ls=":", color="0.5")
plt.ylim(1.e-10*Emax, 1.1*Emax)

plt.savefig("ps.png")


