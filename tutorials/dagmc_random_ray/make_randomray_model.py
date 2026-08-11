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

energy_points = [1.0e-2, 1.0e1]
strengths = [0.0, 0.75]
energy_distribution = openmc.stats.Discrete(x=energy_points, p=strengths)
neutron_source = openmc.IndependentSource(
    energy=energy_distribution,
    space=uniform_dist
)

settings.source = neutron_source

model.settings = settings

# Create voxel plot
plot = openmc.VoxelPlot()
plot.origin = [0, 0, 0]
plot.width = [50, 50, 1]
plot.pixels = [1000, 1000, 1]

model.plots = [plot]

model.export_to_xml()
