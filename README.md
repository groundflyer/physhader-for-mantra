# PhyShader
A set of physical plausible shaders for Mantra renderer.

[![Donate](https://www.paypalobjects.com/webstatic/en_US/btn/btn_donate_74x21.png)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=996RRSDD2C3YQ) [![GitHub release](https://img.shields.io/github/release/groundflyer/physhader-for-mantra.svg)](https://github.com/groundflyer/physhader-for-mantra/releases) [![Houdini Version Compatibilty](https://img.shields.io/badge/houdini-15.5-yellow.svg)](http://www.sidefx.com/index.php?option=com_download&Itemid=208)

## PhySurface VOP
* Energy conserving surface model
* PBR and RayTrace render engines support
* GTR BSDF with anisotropy <sup>[1](#Walter07)</sup> <sup>[2](#Heitz14)</sup> <sup>[3](#Burley12)</sup> _(also avaible as a separate VOP node)_
* Conductor Fresnel
* Volume absorption
* Raytraced subsurface scattering
 * Artist-friendly multiple scattering <sup>[4](#CrBur15)</sup> _(also avaible as a separate VOP node)_
 * Ray-marched single scattering
* Translucency
* Dispersion <sup>[5](#WySlo13)</sup>
* Thin sheet dielectric
* Transparent shadows
* Extra image planes support
 * Per-component image planes
 * Per-light image planes
* Variance anti-aliasing support
* Layered material
* Nesting material <sup>[6](#Schmidt02)</sup>

## PhyVolume VOP
* Color scattering and absorption
* Per-light image planes

##Installation
Copy `vex`,`otls`, `vop` and `gallery` folders into your Houdini home directory`$HOUDINI_USER_PREF_DIR` or `$HIH`.

##Usage
### Quickstart
1. Go to Material Palette
2. Choose PhySurface and move material into the scene
3. Assign material to object
4. Adjust parameters
5. ...
6. RENDER

### Nested dielectrics
[Video Tutorial](https://vimeo.com/180913817)

### Quicktips
1. The preferred Environment Light mode is **Ray tracing background** instead of "Direct Lighting".
2. Indirect lights, such as Ambient Occlusion, Irradiance or Caustic Light do not contribute to subsurface scattering, however Environment Light in Direct Lighting mode does. Try using the combination of two Environment Lights in "Ray Tracing background" and "Direct" modes with object and shader Lightmasks.

## References
<a name="Walter07">1.</a> Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance. Microfacet Models for Refraction through Rough Surfaces. In Proceedings of EGSR 2007.

<a name="Heitz14">2.</a> Eric Heitz, Eugene D'Eon. Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals. Computer Graphics Forum, Wiley-Blackwell, 2014, 33 (4), pp.103-112.

<a name="Burley12">3.</a> Brent Burley. Physically-Based Shading at Disney. 2012.

<a name="CrBur15">4.</a > Per H. Christensen, Brent Burley. An approximate reflectance profile for efficient subsurface scattering. 2015. In ACM SIGGRAPH 2015 Talks (SIGGRAPH '15).

<a name="WySlo13">5.</a> Chris Wyman, Peter-Pike Sloan, and Peter Shirley, Simple Analytic Approximations to the CIE XYZ Color Matching Functions, Journal of Computer Graphics Techniques (JCGT), vol. 2, no. 2, 1-11, 2013.

<a name="Schmidt02">6.</a> Charles M. Schmidt and Brian Budge. Simple Nested Dielectrics In Ray Traced Images. In Journal of Graphics Tools, 7(2), 2002.
