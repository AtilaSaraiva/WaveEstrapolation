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
        i = i + 1
#    start = time()
    gif = manimation.ArtistAnimation(fig, frames, blit=True)
    gif.save(filename, writer='imagemagick', fps=fps)
#    end = time()
#    print('Gif saved; time spent:', end-start, end='.\n\n')
    plt.show()
    plt.close(fig)

    return(P)


def readslice(inputfilename,nx,ny,timeslice):
    f = open(inputfilename,'rb')
    f.seek(8*timeslice*nx*ny)
    field = np.fromfile(f,dtype='float64',count=nx*ny)
    field = np.reshape(field,(nx,ny))
    f.close()
    return field

snaps1 = readslice('snap.ad',200,200,35)
#save_gif(snaps)
