# Definitions

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
The "\f$P\f$ polarisation" refers to the physical set of dipoles parallel to the N-S line; the "\f$Q\f$ polarisation", to the E-W line.

### Local sky coordinates

This is a spherical coordinate system defined with respect to a local observer.
\f$\theta\f$ is the *zenith angle*, i.e. a \textit{colatitude}, with zenith itself therefore defined as \f$\theta = 0\f$ and the horizon as \f$\theta = \pi/2\f$.
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
    \transmat{a}{b}{c}{d} =
    \begin{bmatrix}
        \pd{c}{a} & \pd{c}{b} \\[5 pt]
        \pd{d}{a} & \pd{d}{b}
    \end{bmatrix}
\f]
Among these, the only transformation that is explicitly used in VCSBeam is the transformation between local sky coordinates and celestial sky coordinates, which is a single rotation within the sky plane by the parallactic angle.

