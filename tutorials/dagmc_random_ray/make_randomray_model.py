import openmc

model = openmc.Model()

# materials
# Instantiate some Macroscopic Data
fuel_data = openmc.Macroscopic('fuel')
water_data = openmc.Macroscopic('water')

# Instantiate some Materials and register the appropriate Macroscopic objects
fuel = openmc.Material(name='fuel')
fuel.set_density('macro', 1.0)
fuel.add_macroscopic(fuel_data)
fuel.id = 40

water = openmc.Material(name='water')
water.set_density('macro', 1.0)
water.add_macroscopic(water_data)
water.id = 41

# Instantiate a Materials collection and export to XML
materials = openmc.Materials([fuel, water])
materials.cross_sections = "mgxs.h5"

model.materials = materials

dagmc_univ = openmc.DAGMCUniverse(filename="dagmc.h5m")
geometry = openmc.Geometry(root=dagmc_univ)
model.geometry = geometry

settings = openmc.Settings()
settings.energy_mode = "multi-group"
settings.random_ray['source_shape'] = 'linear'
settings.random_ray['distance_inactive'] = 40.0
settings.random_ray['distance_active'] = 400.0

settings.dagmc = True
settings.batches = 100
settings.inactive = 20
settings.particles = 50000

lower_left  = dagmc_univ.bounding_box.lower_left
upper_right = dagmc_univ.bounding_box.upper_right
uniform_dist = openmc.stats.Box(lower_left, upper_right)
settings.random_ray['ray_source'] = openmc.IndependentSource(space=uniform_dist)

model.settings = settings

model.export_to_xml()
