import openmc

model = openmc.Model()

# materials
u235 = openmc.Material(name="fuel")
u235.add_nuclide('U235', 1.0, 'ao')
u235.set_density('g/cc', 11)
u235.id = 40

water = openmc.Material(name="water")
water.add_nuclide('H1', 2.0, 'ao')
water.add_nuclide('O16', 1.0, 'ao')
water.set_density('g/cc', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')
water.id = 41

mats = openmc.Materials([u235, water])
model.materials = mats

dagmc_univ = openmc.DAGMCUniverse(filename="dagmc.h5m")
geometry = openmc.Geometry(root=dagmc_univ)
model.geometry = geometry

settings = openmc.Settings()
settings.dagmc = True
settings.batches = 100
settings.inactive = 20
settings.particles = 50000

settings.temperature = {'default': 500.0,
                        'method': 'interpolation',
                        'range': (294.0, 3000.0),
                        'tolerance': 1000.0}

settings.source = openmc.IndependentSource(space=openmc.stats.Box([-4., -4., -4.],
                                                                  [ 4.,  4.,  4.]))
model.settings = settings
model.export_to_xml()
model.convert_to_multigroup(method='material_wise',
                            groups='CASMO-2',
                            overwrite_mgxs_library=True,
                            mgxs_path='mgxs.h5',
                            correction = None,
                            source_energy = None,
                            temperatures = None,
                            nparticles = 100_000)

