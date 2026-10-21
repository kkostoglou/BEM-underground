# Chapter 6 wave-propagation validation

The Chapter 6 application parameters were extracted from `Documentation/parts/Chap6_WavePropagation Doc_KostasPhD-2.pdf`.

## Half-space P-wave application

The thesis case (§6.8.1.1) uses a 2.5 m circular cavity whose centre is 5 m below the free surface, E = 20,000 MPa, nu = 0.25, density = 2.55 Mg/m3, a -10 MPa vertical P wave, 20 steps over 0.015 s, and eight quadratic cavity elements.

Full input: [ch6_8_1_1_halfspace_plane_wave.dat](../examples/ch6_8_1_1_halfspace_plane_wave.dat)

A four-element, 20-step validation variant was then completed successfully:

[ch6_8_1_1_halfspace_plane_wave_4elem.dat](../examples/ch6_8_1_1_halfspace_plane_wave_4elem.dat)  
[ch6_8_1_1_halfspace_plane_wave_4elem.npz](../results/ch6_8_1_1_halfspace_plane_wave_4elem.npz)

It produced 8 boundary nodes and 20 time values. Both displacement and traction arrays are finite; maximum absolute displacement is 0.00334194 and maximum absolute stored traction is 10.0 in solver units.

## Explosion application

The thesis explosion application uses the same 2.5 m cavity and material, an underground source at (0, -8) m, equivalent charge weight 188.7 kg, attenuation exponent n = 2.75, coupling factor F = 1.0, and 20 steps over 0.015 s. An explosion-specific Chapter 6 run has not yet been completed because the current wave-propagation path still needs the corresponding free-field/explosion loading and hoop-stress post-processing connected to the graph output.

The existing validated explosion deck `validation/example_decks_1992/original/INFEXP9.DAT` was parsed and executed successfully through the current solver:

[INFEXP9 result](../results/INFEXP9_result.npz)

It produced 16 boundary coordinates and 20 time values, with finite displacement and traction histories. Maximum absolute displacement is 0.00874262 and maximum absolute stored traction is 37.7950 in solver units. This validates the explosion loading path numerically, although the historical thesis hoop-stress curves still require stress reconstruction for direct graph comparison.

## Validation status

The four-element, 20-step half-space case validates that the current time-domain input path accepts the Chapter 6 geometry and plane-wave loading. Thesis-level validation remains pending until the full run completes and boundary tangential stress is reconstructed for comparison with Figures 6.18 and 6.19.
