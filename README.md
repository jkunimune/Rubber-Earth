# The Danseiji projection

The Dansēji (or Danseiji; /dɑn.ˈseɪ.dʒi/; literally "Elastic Earth") projections are the optimal map projections,
in the most literal sense of the word. Each one is represented by a fine mesh of vertices
that represent an elastic spherical surface. That surface is allowed to come to rest, thus
reducing the distortions of size and shape throughout the mesh.

For a more complete explanation, see
[this blog post](https://kunimune.home.blog/2019/11/07/introducing-the-danseiji-projections/).
For the newer map projections based on the same principle, see my repository for
[the Elastic projections](https://github.com/jkunimune/elastic).

## Use

All of the source code and output data files are released to the public domain (see [LICENSE](https://github.com/jkunimune/Rubber-Earth/blob/master/LICENSE) for more information). The Dansēji projections are defined by the CSV files found in [output/](https://github.com/jkunimune15/Rubber-Earth/tree/master/output). An example Python script using these to create SVG maps can be found in [src/example/](https://github.com/jkunimune15/Rubber-Earth/tree/master/src/example).

## Theory

Each Dansēji projection is created, in theory, with an inifintely thin spherical shell of
an elastic material, representing the Earth's surface, which is then constrained into a
plane. The initial configuration is either an equatorial Hammer projection, for a
more conventional map, or an oblique azimuthal equidistant projection, for a more
optimal one. This surface is then deformed in order to mimimize the elastic
potential energy under a Neo-Hookean model. In practice, I do this using finite element
analysis and an L-BFGS scheme.

## Credits

* [The GeoTools library](http://docs.geotools.org/)
* [The Natural Earth Data](https://www.naturalearthdata.com/)
* [The NASA Earth Observatory](https://neo.sci.gsfc.nasa.gov/)
* [NASA Socioeconomic Data and Applications Center](https://sedac.ciesin.columbia.edu/)
