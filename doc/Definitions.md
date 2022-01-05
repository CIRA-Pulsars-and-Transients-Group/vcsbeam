# Definitions {#definitions}

[TOC]

## Coordinate systems

There are three coordinate systems in use throughout VCSBeam:
 - Instrumental (\f$p\f$,\f$q\f$)
 - Local sky coordinates (\f$\theta\f$,\f$\phi\f$)
 - Celestial sky coordinates (\f$x\f$,\f$y\f$)
They are illustrated in the following figure:

\image html coords.png width=800
\image latex coords.png width=\textwidth

### Instrumental coordinates

This is a Cartesian coordinate system aligned with local (ground) compass directions.
Positive \f$p\f$ points towards local North, and positive \f$q\f$ towards local East.
The \f$P\f$ polarisation refers to the physical set of dipoles parallel to the N-S line; the \f$Q\f$ polarisation, to the E-W line.

### Local sky coordinates

This is a spherical coordinate system defined with respect to a local observer.
\f$\theta\f$ is the *zenith angle*, i.e. a *colatitude*, with zenith itself therefore defined as \f$\theta = 0\f$ and the horizon as \f$\theta = \pi/2\f$.
\f$\phi\f$ is the *azimuth*, and we define \f$\phi = 0\f$ in the North direction, with positive azimuth moving clockwise as viewed from above (i.e. N&rarr;E&rarr;S&rarr;W&rarr;N).
Moreover, the *elevation* is denoted by the symbol \f$\tilde{\theta}\f$, and is related to the zenith angle by
\f[
    \tilde{\theta} = \frac{\pi}{2} - \theta.
\f]

### Celestial sky coordinates

This is a spherical coordinate system defined with respect to the celestial sphere.
\f$x\f$ is the *declination* (Dec) and \f$y\f$ is the *right ascension* (RA).

## Coordinate transformations

All coordinate transformations can be effected by applying the appropriate Jacobian matrix for the desired transformation.
In this documentation, a boldface \f${\bf P}\f$ will always be used to denote transformation matrices.
For general coordinates \f$(a,b)\f$ and \f$(c,d)\f$,
\f[
    {\bf P}_{(a,b)\rightarrow(c,d)} \equiv
    \begin{bmatrix}
        P_{ca} & P_{cb} \\
        P_{da} & P_{db}
    \end{bmatrix}
    \equiv
    \begin{bmatrix}
        \frac{\partial c}{\partial a} & \frac{\partial c}{\partial b} \\[5 pt]
        \frac{\partial d}{\partial a} & \frac{\partial d}{\partial b}
    \end{bmatrix}
\f]
Among these, the only transformation that is explicitly used in VCSBeam is the transformation between local sky coordinates and celestial sky coordinates, which is a single rotation within the sky plane by the parallactic angle.

### Parallactic angle correction {#parallacticangle}

The parallactic angle correction is a transformation between local sky coordinates and celestial sky coordinates.
The parallactic angle itself, \f$\chi\f$, is defined as the position angle of local zenith with respect to the North Celestial Pole as subtended at a given source, illustrated in the following figure.

\image html skyangles.png width=300
\image latex skyangles.png width=0.4\textwidth

The transformation \f$(x,y)\rightarrow(\theta,\phi)\f$ is therefore a counterclockwise rotation by \f$\pi - \chi\f$.
(A counterclockwise rotation of a given vector is equivalent to a clockwise rotation of the coordinate axes.)
This is the rotation
\f[
    {\bf P}_\text{pa}
        = \begin{bmatrix}
            P_{\theta x} & P_{\theta y} \\
            P_{\phi x}   & P_{\phi y}
        \end{bmatrix}
        = \begin{bmatrix}
            \cos\left(\pi - \chi\right) & -\sin\left(\pi - \chi\right) \\
            \sin\left(\pi - \chi\right) &  \cos\left(\pi - \chi\right)
        \end{bmatrix}
        = \begin{bmatrix}
            -\cos\chi & -\sin\chi \\
             \sin\chi & -\cos\chi
        \end{bmatrix}.
\f]

In VCSBeam, the parallactic angle is calculated (via the function `palPa`) in the [Starlink/pal library](https://github.com/Starlink/pal) by the spherical triangle identity
\f[
    \tan \chi = \frac{\cos \lambda \sin H}{\sin\lambda \cos x - \cos \lambda \sin x \cos H},
\f]
where \f$\lambda\f$ is the latitude of the observer, \f$H\f$ is the hour angle of the source, and \f$x\f$ is the declination.
The latitude of the MWA is \f$\lambda_\text{MWA} = -0.4660608448386394\f$ rad, defined in the [mwalib](https://github.com/MWATelescope/mwalib) library.

## Data dimensions {#datadimensions}

All data are functions of time (*t*), frequency (*f*), antenna (*a*), and/or polarisation (*p*).
Accordingly, the size of data arrays are expressed in terms of \f$N_t\f$, \f$N_f\f$, \f$N_a\f$, and \f$N_p\f$.
The order of the dimensions indicates the layout of the data arrays in memory, with the first dimension being the slowest changing dimension, and so on until the last dimension, which is the fastest changing.
For example, the instrumental calibration matrices (\f${\bf D}\f$) have dimensions
\f[
N_a \times N_f \times N_p \times N_p,
\f]
and a single element is specified with the index notation, \f$D_{a,f,p_1,p_2}\f$.
*Jones vectors*, indicated with a lower case Latin letter in **bold**, include a polarisation dimension implicitly, e.g.
\f[ {\bf e} = \begin{bmatrix} e_x \\ e_y \end{bmatrix}, \f]
where \f$x\f$ and \f$y\f$ are examples of a specific polarisation basis.
Similary, *Jones matrices*, represented with **bold** uppercase Latin letters, have two implied polarisation dimensions, e.g.
\f[ {\bf J} = \begin{bmatrix} J_{xx} & J_{xy} \\ J_{yx} & J_{yy} \end{bmatrix}. \f]

The exception to the above system is the data layout of the legacy recombined format, in which the combinations of antennas and polarisations are ordered in a "mixed" way, so that it cannot be said that either antennas or polarisations change faster than the other.
For this case, we describe a particular antenna-polarisation combination as an "input", and denote it with the subscript \f$i\f$, with \f$N_i = N_a \times N_p\f$ (without the implied hierarchical ordering).

Full Stokes parameters are the set \f$\{I, Q, U, V\}\f$.
Accordingly, an array containing these "post-detection" products may include the dimension \f$N_s\f$, which always equals exactly 4.

Each MWA tile contains \f$N_d = 16\f$ dipoles.

Multiple pointings can also be processed simultaneously, and some output arrays therefore include dimension \f$N_b\f$.

The complete list, in alphabetical order of subscripts, is
| Symbol    | Description                               | Typical values |
| :-------: | :---------------------------------------- | :------------: |
| \f$N_a\f$ | Number of tiles/antennas                  | 128            |
| \f$N_b\f$ | Number of beams/pointings                 | 1              |
| \f$N_d\f$ | Number of dipoles per tile                | 16             |
| \f$N_f\f$ | Number of frequencies per coarse channel  | 128            |
| \f$N_p\f$ | Number of polarisations                   | 2              |
| \f$N_s\f$ | Number of Stokes parameters               | 4              |
| \f$N_t\f$ | Number of timesteps/samples per second    | 10000          |

## Comparison of notation in other documents

### Coordinate systems

| Document           |         |         |              |            |         |         |
| :----------------- | :-----: | :-----: | :----------: | :--------: | :-----: | :-----: |
| *This document*    | \f$p\f$ | \f$q\f$ | \f$\theta\f$ | \f$\phi\f$ | \f$x\f$ | \f$y\f$ |
| MWA metafits files | Y       | X       | -            | -          | -       | -       |

### Jones matrices

[Sokolowski2017]: https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/calibration-and-stokes-imaging-with-full-embedded-element-primary-beam-model-for-the-murchison-widefield-array/FBA84B9EB94000BD6258A8F75840C476
[Ord2019]: https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/abs/mwa-tiedarray-processing-i-calibration-and-beamformation/E9A7A9981AE9A935C9E08500CA6A1C1E

| Document                                    |               |               |               |
| :------------------------------------------ | :-----------: | :-----------: | :-----------: |
| *This document*                             | \f${\bf J}\f$ | \f${\bf D}\f$ | \f${\bf B}\f$ |
| [Sokolowski et al., (2017)][Sokolowski2017] | \f${\bf J}\f$ | \f${\bf G}\f$ | \f${\bf E}\f$ |
| [Ord et al., (2019)][Ord2019]               | \f${\bf J}\f$ | \f${\bf D}\f$ | \f${\bf B}\f$ |
