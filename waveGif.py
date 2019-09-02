import wave

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import seaborn as sns
from numba import jit, prange, typeof
from os import system

def save_gif(P, name='animation',
    fps=16,
    vmin=None, vmax=None,
    vwin=None,
    cmap='rainbow', dpi=150):
    '''
    Saves a [t,z,x] array as a gif.
    '''
    filename = name + '.gif'

    frames_div = 5

    nx, nz, nt = P.shape
    fig = plt.figure(dpi=dpi)
    fig.set_size_inches(2*nx/np.linalg.norm([nx,nz]),
               2*nz/np.linalg.norm([nx,nz]))

    ax = fig.add_axes([0, 0, 1, 1], frameon=False, aspect=1)
    ax.set_xticks([])
    ax.set_yticks([])

    if vwin:
        for coord in ['xi', 'zi', 'xf', 'zf']:
            if not coord in vwin.keys(): vwin[coord] = None
        vview = P[:,vwin['zi']:vwin['zf'],vwin['xi']:vwin['xf']]
        vmin = vview.min()
        vmax = vview.max()
    else:
        if not vmin: vmin = P.min()
        if not vmax: vmax = P.max()

    frames = []
    i=0
    while i<nt:
        frames.append([ax.imshow(P[:,:,i], cmap, animated=True,vmin=vmin,vmax=vmax,
                interpolation='nearest')])
        i = i + frames_div
#    start = time()
    gif = manimation.ArtistAnimation(fig, frames, blit=True)
    gif.save(filename, writer='imagemagick', fps=fps)
#    end = time()
#    print('Gif saved; time spent:', end-start, end='.\n\n')
    plt.show()
    plt.close(fig)

    return(P)

L  = 10
Nx = 200
Nz = 200
Nt = 300
h1 = 1000
h2 = 1000
dx = 3.0
dz = 3.0
#dt = (4000*(dx**(-2)+dz**(-2)))**(-0.5)
dt = 0.02

wavelet = np.empty(10)
for i in range(len(wavelet)):
    #wavelet[i] = 20.0*np.sin(20.0*(i-1)*dt)
    wavelet[i] = 1 + (i+1)**(-1) - (i+1)**(-2) + (i+1)**(-3)
#h1A = int(h1/dz)
#h2A = int(h2/dz)
h1A = int(Nz/3)
h2A = int(Nz/3)

campoVel                = np.zeros((Nz,Nx))
campoVel[0:h1A,:]       = 2000.0
campoVel[h1A:h2A+h1A,:] = 3500.0
campoVel[h2A+h1A:Nz,:]  = 4000.0

fontes = np.empty((2,1))
fontes[:,0] = np.array([1,int(Nx/2)])

snapCube = wave.diffinitas.waveestrap (Nt,dx,dz,dt,wavelet,campoVel,fontes,nx=Nx,nz=Nz)

"""
fig, ax = plt.subplots()
cmap='rainbow'
parada = ax.imshow(snapCube[:,:,50], cmap, animated=True,
            interpolation='nearest')
frames = []
frames_div = 5
for frame in snapCube[1::frames_div]:
    frames.append(
    [ax.imshow(frame, cmap, animated=True,
            interpolation='nearest')])

frame = []
for i in range(Nt):
    frames.append([ax.imshow(snapCube[:,:,i], cmap, animated=True,
            interpolation='nearest')])
gif = manimation.ArtistAnimation(fig, frames, blit=True)

plt.show()

"""
save_gif(snapCube,name="wave",fps=24)
