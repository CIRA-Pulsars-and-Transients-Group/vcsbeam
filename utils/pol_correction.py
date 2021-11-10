#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

def phase_correction(f_MHz, slope, offset):
    return slope*f_MHz*1e6 + offset


def load_pdv(filename):
    pdv = np.loadtxt('vela.pdv')

    ntime = int(pdv[-1, 0] + 1)
    nchan = int(pdv[-1, 1] + 1)
    nbins = int(pdv[-1, 2] + 1)

    S = {"I": np.reshape(pdv[:, 3], (ntime, nchan, nbins)),
         "Q": np.reshape(pdv[:, 4], (ntime, nchan, nbins)),
         "U": np.reshape(pdv[:, 5], (ntime, nchan, nbins)),
         "V": np.reshape(pdv[:, 6], (ntime, nchan, nbins))}

    return S, ntime, nchan, nbins

def load_beam(filename, ntime, nchan):
    beam = np.loadtxt('beam.txt')

    J = np.array([[[beam[i,10] + beam[i,11]*1j, beam[i,12] + beam[i,13]*1j], [beam[i,14] + beam[i,15]*1j, beam[i,16] + beam[i,17]*1j]] for i in range(beam.shape[0])])
    J = J.reshape((nchan, ntime, 2, 2))

    freqs_MHz = beam[::ntime,1]
    return J, freqs_MHz

class InteractiveSlopePlot:

    def __init__(self, S, Qs, freqs_MHz, slope, offset):

        self.S = S
        self.Qs = Qs
        self.active = False
        self.points_plot = None
        self.lines_plot = None
        self.slope = slope
        self.offset = offset
        self.freqs_MHz = freqs_MHz
        self.click_number = 0

        self.fig = plt.figure(constrained_layout=True)
        gs = self.fig.add_gridspec(2, 6)
        ax00 = self.fig.add_subplot(gs[0, 0])
        ax01 = self.fig.add_subplot(gs[0, 1])
        ax02 = self.fig.add_subplot(gs[0, 2])
        ax03 = self.fig.add_subplot(gs[0, 3])
        ax10 = self.fig.add_subplot(gs[1, 0])
        ax11 = self.fig.add_subplot(gs[1, 1])
        ax12 = self.fig.add_subplot(gs[1, 2])
        ax13 = self.fig.add_subplot(gs[1, 3])
        self.ax = self.fig.add_subplot(gs[:, 4:])

        self.axs_orig = [ax00, ax01, ax02, ax03]
        self.axs_corr = [ax10, ax11, ax12, ax13]

        self.apply_correction()
        self.plot_phase_correction()

        self.plot_all_stokes("orig")
        self.plot_all_stokes("corr")

        # Turn interactivity on
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_button_press_event)
        self.cid = self.fig.canvas.mpl_connect('axes_enter_event', self.on_axes_enter_event)
        self.cid = self.fig.canvas.mpl_connect('axes_leave_event', self.on_axes_leave_event)

        plt.show()

    def apply_correction(self):

        # Shorthand variables
        S         = self.S
        freqs_MHz = self.freqs_MHz
        Qs        = self.Qs
        slope     = self.slope
        offset    = self.offset

        Zs = np.array([[[1, 0], [0, np.exp(1j*ph)]] for ph in phase_correction(freqs_MHz, slope, offset)])
        Zinvs = np.array([np.linalg.pinv(Zs[i,:,:]) for i in range(Zs.shape[0])])
        Qinvs = np.array([[np.linalg.pinv(Qs[i,j,:,:]) for j in range(Qs.shape[1])] for i in range(Qs.shape[0])])

        QZQs = np.array([[Qs[i,j,:,:] @ Zinvs[i,:,:] @ Qinvs[i,j,:,:] for j in range(Qs.shape[1])] for i in range(Qs.shape[0])])

        # Now convert these to Mueller matrices
        # (see https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470060193.app4)
        # ... but I'm going to do it long-hand
        Lxx = QZQs[:,:,0,0]
        Lxy = QZQs[:,:,0,1]
        Lyx = QZQs[:,:,1,0]
        Lyy = QZQs[:,:,1,1]

        MII = Lxx*np.conj(Lxx) + Lxy*np.conj(Lxy) + Lyx*np.conj(Lyx) + Lyy*np.conj(Lyy)
        MIQ = Lxx*np.conj(Lxx) - Lxy*np.conj(Lxy) + Lyx*np.conj(Lyx) - Lyy*np.conj(Lyy)
        MIU =  2*np.real(Lxx*np.conj(Lxy) + Lyx*np.conj(Lyy))
        MIV = -2*np.imag(Lxx*np.conj(Lxy) + Lyx*np.conj(Lyy))

        MQI = Lxx*np.conj(Lxx) + Lxy*np.conj(Lxy) - Lyx*np.conj(Lyx) - Lyy*np.conj(Lyy)
        MQQ = Lxx*np.conj(Lxx) - Lxy*np.conj(Lxy) - Lyx*np.conj(Lyx) + Lyy*np.conj(Lyy)
        MQU =  2*np.real(Lxx*np.conj(Lxy) - Lyx*np.conj(Lyy))
        MQV = -2*np.imag(Lxx*np.conj(Lxy) - Lyx*np.conj(Lyy))

        MUI =  2*np.real(Lxx*np.conj(Lyx) + Lxy*np.conj(Lyy))
        MUQ =  2*np.real(Lxx*np.conj(Lyx) - Lxy*np.conj(Lyy))
        MUU =  2*np.real(Lxx*np.conj(Lyy) + Lxy*np.conj(Lyx))
        MUV = -2*np.imag(Lxx*np.conj(Lyy) - Lxy*np.conj(Lyx))

        MVI =  2*np.imag(Lxx*np.conj(Lyx) + Lxy*np.conj(Lyy))
        MVQ =  2*np.imag(Lxx*np.conj(Lyx) - Lxy*np.conj(Lyy))
        MVU =  2*np.imag(Lxx*np.conj(Lyy) + Lxy*np.conj(Lyx))
        MVV =  2*np.real(Lxx*np.conj(Lyy) - Lxy*np.conj(Lyx))

        M = np.array([[[[MII[i,j], MIQ[i,j], MIU[i,j], MIV[i,j]],
                        [MQI[i,j], MQQ[i,j], MQU[i,j], MQV[i,j]],
                        [MUI[i,j], MUQ[i,j], MUU[i,j], MUV[i,j]],
                        [MVI[i,j], MVQ[i,j], MVU[i,j], MVV[i,j]]] for j in range(MII.shape[1])] for i in range(MII.shape[0])])

        # With this Mueller matrix, now apply it!
        Svecs = np.array([[[[S["I"][i,j,k], S["Q"][i,j,k], S["U"][i,j,k], S["V"][i,j,k]] for i in range(S["I"].shape[0])] for j in range(S["I"].shape[1])] for k in range(S["I"].shape[2])])

        Sprime = np.array([[[M[i,j,:,:] @ Svecs[k,i,j,:] for k in range(Svecs.shape[0])] for i in range(Svecs.shape[1])] for j in range(Svecs.shape[2])])

        self.Scorrected = {"I": np.real(Sprime[:,:,:,0]),
                           "Q": np.real(Sprime[:,:,:,1]),
                           "U": np.real(Sprime[:,:,:,2]),
                           "V": np.real(Sprime[:,:,:,3])}

    def on_button_press_event(self, event):
        if self.active == True:
            if self.click_number == 0:

                self.point1 = np.array([event.xdata*1e6, event.ydata])

                if self.points_plot is None:
                    self.points_plot, = self.ax.plot([self.point1[0]/1e6], [self.point1[1]], 'rx')
                else:
                    self.points_plot.set_data([self.point1[0]/1e6], [self.point1[1]])

                self.fig.canvas.draw()

            elif self.click_number == 1:

                self.point2 = np.array([event.xdata*1e6, event.ydata])
                diffs = self.point2 - self.point1
                self.slope = np.deg2rad(diffs[1])/diffs[0]
                self.offset = np.deg2rad(self.point1[1]) - self.slope*self.point1[0]

                # Recalculate the corrected Stokes
                self.apply_correction()
                self.plot_all_stokes("corr")
                self.plot_phase_correction()
                self.points_plot.set_data([], [])

                self.fig.canvas.draw()

            self.click_number = (self.click_number + 1) % 2

    def on_axes_enter_event(self, event):
        if event.inaxes == self.ax:
            self.active = True

    def on_axes_leave_event(self, event):
        if event.inaxes == self.ax:
            self.active = False

    def plot_phase_correction(self):
        ph = phase_correction(self.freqs_MHz, self.slope, self.offset)
        ph = (ph + np.pi) % (2 * np.pi) - np.pi
        if self.lines_plot is None:
            self.lines_plot, = self.ax.plot(self.freqs_MHz, np.rad2deg(ph))
        else:
            self.lines_plot.set_data(self.freqs_MHz, np.rad2deg(ph))
        self.ax.set_xlabel("Frequency (MHz)")
        self.ax.set_ylabel("Phase (deg)")
        self.ax.set_ylim(-180, 180)
        self.ax.set_title("Slope = {} rad/Hz\nOffset = {}".format(self.slope, self.offset))

    def plot_stokes(self, s, ax, t=None):
        if t is None:
            image = np.sum(s, axis=0)
        else:
            image = s[t,:,:]

        if self.freqs_MHz is None:
            extent=None
        else:
            df = freqs_MHz[1] - freqs_MHz[0]
            extent=(-0.5, image.shape[-1] - 0.5, self.freqs_MHz[0] - df/2, self.freqs_MHz[-1]  + df/2)
        ax.imshow(image, origin='lower', interpolation='none', aspect='auto', cmap='hot', extent=extent)

    def plot_all_stokes(self, orig_or_corr="orig"):
        if orig_or_corr == "orig":
            S = self.S
            axs = self.axs_orig
        else:
            S = self.Scorrected
            axs = self.axs_corr

        ss = list(S.keys())
        for i in range(len(ss)):
            s  = S[ss[i]]
            ax = axs[i]
            self.plot_stokes(s, ax)

            # Set up axes
            if i > 0:
                ax.set_yticks([])
            else:
                if self.freqs_MHz is None:
                    ax.set_ylabel("Frequency channel number")
                else:
                    ax.set_ylabel("Frequency (MHz)")

            ax.set_xlabel("Phase bin")

if __name__ == '__main__':
    pdvfile = sys.argv[1]
    S, ntime, nchan, nbins = load_pdv(pdvfile)

    beamfile = sys.argv[2]
    Qs, freqs_MHz = load_beam(beamfile, ntime, nchan)

    if len(sys.argv) > 3:
        slope = float(sys.argv[3])
    else:
        slope = 0.0

    if len(sys.argv) > 4:
        offset = float(sys.argv[4])
    else:
        offset = 0.0

    slope_plot = InteractiveSlopePlot(S, Qs, freqs_MHz, slope, offset)
