# The Dansēji ("Elastic Earth") projection

For centuries, cartographers have wrestled with the challenges of making accurate maps.
They were limited by not only the mathematical impossibility of accurately portraying a
curved surface on a plane, but also by their own computational limitations. In Mercator's
day, creating a single conformal cylindrical map required performing a Riemann sum for
every single latitude. Even today, some of the most mathematically complex map
projections cling to simple geometric shapes, like the straight parallels of Tobler or
the tetrahedrons of AuthaGraph. Modern cartographers work to minimise distortion and
constrain it to the least important regions of the globe -- usually the oceans -- within
the constraints of these basic formulae through a combination of trial-and-error,
combining and varying existing projections, and optimisation by advanced calculus.

And yet, this is the twenty-first century. The Mayan apocalypse has come and gone. The
singularity is nigh on the horizon. Cartographers need not bow to these computational
constraints any longer. The time has come to shake free of these arbitrary lines of
latitude and longitude, to throw off these circles and hyperellipses and power
functions, to deny control over our maps to these unchallenged assumptions that physical
area is the most important quantity of a country, that Antarctica is unimportant and
negligible, that the International Dateline is the end of the world! The time has come
for the Dansēji projections.

The Dansēji projections are the optimal map projections, in the most literal sense of
the word. Each one is represented by a fine mesh of vertices that represent an elastic
spherical surface. That surface is allowed to come to rest, thus reducing the
distortions of size and shape throughout the mesh. The core of the Dansēji projection's
strength comes from the fact that it _has no equations_. It is free to be anything or
anyone, regardless of how ungeometrical or irregular that may be. This family of maps
can be easily adapted and applied to any function. Comparing areas? Showing statistics?
Teaching children? There's a Dansēji projection for you!

Unlike the Cahill, Dymaxion, and AuthaGraph projections, whose equations are complex
and poorly-documented (or even obscured from the public under copyright in AuthaGraph's
case), the Dansēji projection has its coordinates published in this repository under
(I'll do this later) in a simple format, for anyone with basic computing knowledge and
an understanding of linear interpolation to use for their own renditions.

## Theory

The Dansēji projection is created, in theory, with an inifintely thin spherical shell of
an elastic material, representing the Earth's surface, which is then constrained into a
plane. The initial configuration is either an equatorial sinusoidal projection, for a
more conventional map, or an oblique azimuthal equidistant projection, for a more
optimal one. This mesh is then deformed in order to mimimise the elastic
potential energy under a Neo-Hookean model. In practise, I do this using finite element
analysis and L-BFGS approximation.

## Credits

* [The GeoTools library](http://docs.geotools.org/)
* [The Natural Earth Data](https://www.naturalearthdata.com/)
* [The NASA Earth Observatory](https://neo.sci.gsfc.nasa.gov/)
